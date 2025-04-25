#!/bin/bash

echo "Script iniciado."
echo "Se você colou este script no terminal, pressione Enter para continuar (isso evita execução automática acidental)."
read

# Informações e verificação de imagens Docker
echo ""
echo "Este script executará análises de qualidade utilizando Docker."
echo "As imagens necessárias são:"
echo "- biocontainers/fastqc:v0.11.9_cv8"
echo "- ewels/multiqc"

MISSING_IMAGES=()

if ! docker image inspect biocontainers/fastqc:v0.11.9_cv8 &>/dev/null; then
  MISSING_IMAGES+=("biocontainers/fastqc:v0.11.9_cv8")
fi

if ! docker image inspect ewels/multiqc &>/dev/null; then
  MISSING_IMAGES+=("ewels/multiqc")
fi

if [[ ${#MISSING_IMAGES[@]} -gt 0 ]]; then
  echo ""
  echo "As seguintes imagens Docker estão faltando:"
  for img in "${MISSING_IMAGES[@]}"; do echo "- $img"; done
  echo ""
  echo "Por favor, baixe as imagens necessárias com:"
  for img in "${MISSING_IMAGES[@]}"; do echo "  docker pull $img"; done
  echo ""
  echo "Após instalar as imagens, execute o script novamente."
  exec bash
else
  echo "Verificação concluída: todas as imagens Docker necessárias estão disponíveis."
fi

# Pergunta sobre uso de --user
echo ""
echo "Deseja rodar o Docker com um usuário específico (--user UID)?"
echo "Digite o UID desejado (ex: 1006), ou 'no' para continuar normalmente:"
read CUSTOM_UID
if [[ "$CUSTOM_UID" == "exit" ]]; then exec bash; fi
if [[ "$CUSTOM_UID" =~ ^[0-9]+$ ]]; then
  USE_USER="--user $CUSTOM_UID"
  echo "Docker será executado com: $USE_USER"
elif [[ "$CUSTOM_UID" == "no" ]]; then
  USE_USER=""
  echo "Docker será executado com usuário padrão (root dentro do container)."
else
  echo "Entrada inválida. Digite um número de UID ou 'no'."
  exec bash
fi

echo ""
echo "Informe o caminho do diretório com os arquivos *_R1_001.fastq.gz* e *_R2_001.fastq.gz* (ou digite 'exit' para sair):"
read INPUT_DIR
if [[ "$INPUT_DIR" == "exit" ]]; then exec bash; fi
if [[ ! -d "$INPUT_DIR" ]]; then
  echo "Erro: diretório '$INPUT_DIR' não encontrado."
  exec bash
fi

echo "Listando arquivos *_val_1.fq.gz e *_val_2.fq.gz em '$INPUT_DIR'..."
FILES=($(find "$INPUT_DIR" -maxdepth 1 \( -name "*_val_1.fq.gz" -o -name "*_val_2.fq.gz" \) -printf "%f\n"))
for file in "${FILES[@]}"; do echo "$file"; done

echo "Deseja processar todos os arquivos listados acima? (yes/no ou 'exit')"
read PROCESS_ALL
if [[ "$PROCESS_ALL" == "exit" ]]; then exec bash; fi

if [[ "$PROCESS_ALL" != "yes" ]]; then
  echo "Insira os nomes EXATOS dos arquivos que deseja processar, separados por espaço (ou 'exit' para sair):"
  read -a FILES
  if [[ "${FILES[0]}" == "exit" ]]; then exec bash; fi
fi

# Criar diretório de saída
TIMESTAMP=$(date +"%d-%m-%Y_%Hh%Mm")
OUTPUT_DIR="$INPUT_DIR/output_qc_$TIMESTAMP"
mkdir -p "$OUTPUT_DIR"

echo "Os arquivos de saída serão salvos em: $OUTPUT_DIR"
echo "Deseja continuar? (yes/no ou 'exit')"
read CONTINUE
if [[ "$CONTINUE" == "exit" ]]; then exec bash; fi
if [[ "$CONTINUE" != "yes" ]]; then
  echo "Processo cancelado."
  exec bash
fi

# Associar arquivos por paciente
declare -A PATIENTS

for FILE in "${FILES[@]}"; do
  if [[ "$FILE" == *_val_1.fq.gz ]]; then
    ID="${FILE%_val_1.fq.gz}"
    PATIENTS["$ID,R1"]="$FILE"
  elif [[ "$FILE" == *_val_2.fq.gz ]]; then
    ID="${FILE%_val_2.fq.gz}"
    PATIENTS["$ID,R2"]="$FILE"
  fi
done

echo ""
echo "Iniciando análise de qualidade..."

for ID in "${!PATIENTS[@]}"; do
  [[ "$ID" =~ ,R1$ ]] || continue
  BASE_ID="${ID%,R1}"
  R1_FILE="${PATIENTS[$BASE_ID,R1]}"
  R2_FILE="${PATIENTS[$BASE_ID,R2]}"

  if [[ -z "$R1_FILE" || -z "$R2_FILE" ]]; then
    echo "Aviso: par incompleto para '$BASE_ID'. Pulando."
    continue
  fi

  # Verificação de nomes invertidos
  if [[ "$R1_FILE" == *R2* || "$R2_FILE" == *R1* ]]; then
    echo ""
    echo "Atenção: os nomes dos arquivos indicam que R1 e R2 podem estar invertidos."
    echo "Detectado:"
    echo "  R1_FILE: $R1_FILE"
    echo "  R2_FILE: $R2_FILE"
    echo "Deseja continuar mesmo assim? (yes/no)"
    read CONFIRM_INVERTED
    if [[ "$CONFIRM_INVERTED" != "yes" ]]; then
      echo "Processo interrompido para o par $BASE_ID."
      continue
    fi
  fi

  echo ""
  echo "Processando paciente: $BASE_ID"
  echo "Arquivo R1: $R1_FILE"
  echo "Arquivo R2: $R2_FILE"

  # FastQC para R1
  docker run --rm \
    -v "$INPUT_DIR":/data \
    -v "$OUTPUT_DIR":/output \
    $USE_USER \
    biocontainers/fastqc:v0.11.9_cv8 \
    fastqc /data/"$R1_FILE" -o /output

  mv "$OUTPUT_DIR/${R1_FILE%.fq.gz}_fastqc.html" "$OUTPUT_DIR/${BASE_ID}_R1_fastqc.html" 2>/dev/null
  mv "$OUTPUT_DIR/${R1_FILE%.fq.gz}_fastqc.zip" "$OUTPUT_DIR/${BASE_ID}_R1_fastqc.zip" 2>/dev/null

  # FastQC para R2
  docker run --rm \
    -v "$INPUT_DIR":/data \
    -v "$OUTPUT_DIR":/output \
    $USE_USER \
    biocontainers/fastqc:v0.11.9_cv8 \
    fastqc /data/"$R2_FILE" -o /output

  mv "$OUTPUT_DIR/${R2_FILE%.fq.gz}_fastqc.html" "$OUTPUT_DIR/${BASE_ID}_R2_fastqc.html" 2>/dev/null
  mv "$OUTPUT_DIR/${R2_FILE%.fq.gz}_fastqc.zip" "$OUTPUT_DIR/${BASE_ID}_R2_fastqc.zip" 2>/dev/null

  # MultiQC (analisando R1 e R2 juntos)
  TEMP_MULTIQC_DIR="$OUTPUT_DIR/multiqc_temp_$BASE_ID"
  mkdir -p "$TEMP_MULTIQC_DIR"
  cp "$OUTPUT_DIR/${BASE_ID}_R1_fastqc.zip" "$TEMP_MULTIQC_DIR"
  cp "$OUTPUT_DIR/${BASE_ID}_R2_fastqc.zip" "$TEMP_MULTIQC_DIR"

  docker run --rm \
    -v "$TEMP_MULTIQC_DIR":/data \
    -v "$OUTPUT_DIR":/output \
    $USE_USER \
    ewels/multiqc \
    multiqc /data --outdir /output --filename "${BASE_ID}_multiqc_report.html"

  rm -r "$TEMP_MULTIQC_DIR"
done

echo ""
echo "Análise de qualidade concluída."
echo "Todos os arquivos foram salvos em: $OUTPUT_DIR"
exec bash
