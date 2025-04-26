#!/bin/bash

# Início
echo "========================================="
echo "        Script para rodar MultiQC         "
echo "========================================="
echo ""
echo "Se colou esse script no terminal, pressione Enter para continuar."
read

# Pergunta se quer usar UID
echo "Deseja rodar o Docker com um usuário específico (--user UID)?"
echo "Digite o número do UID (ex: 1006), ou 'no' para seguir como root, ou 'exit' para sair."
read CUSTOM_UID
if [[ "$CUSTOM_UID" == "exit" ]]; then exit 0; fi

if [[ "$CUSTOM_UID" =~ ^[0-9]+$ ]]; then
  USE_USER="--user $CUSTOM_UID"
  echo "Docker será executado com: $USE_USER"
elif [[ "$CUSTOM_UID" == "no" ]]; then
  USE_USER=""
  echo "Docker será executado como root dentro do container."
else
  echo "Entrada inválida. Encerrando."
  exit 1
fi

# Onde estão os arquivos
echo "Agora informe o caminho COMPLETO onde estão os arquivos (todos devem estar no mesmo diretório)."
echo "Para garantir, acesse a pasta onde estão os arquivos, digite 'pwd', copie o caminho e cole aqui."
read INPUT_DIR
if [[ "$INPUT_DIR" == "exit" ]]; then exit 0; fi

if [[ ! -d "$INPUT_DIR" ]]; then
  echo "Erro: Diretório '$INPUT_DIR' não encontrado."
  exit 1
fi

# Mostrar os arquivos
echo "Listando todos os arquivos em '$INPUT_DIR'..."
FILES_AVAILABLE=($(ls "$INPUT_DIR"))
for file in "${FILES_AVAILABLE[@]}"; do
  echo "- $file"
done

# Receber nomes dos arquivos para MultiQC
echo "Agora copie e cole os nomes dos arquivos que deseja incluir no MultiQC, separados por espaço."
echo "ATENÇÃO: Informe apenas os arquivos .zip (exemplo abaixo):"
echo "paciente1_R1_fastqc.zip paciente1_R2_fastqc.zip paciente2_R1_fastqc.zip paciente2_R2_fastqc.zip"
echo "(Use exatamente os nomes mostrados acima, respeitando maiúsculas, minúsculas e extensões!)"
read -a FILES_SELECTED

# Validação básica
if [[ -z "${FILES_SELECTED[*]}" ]]; then
  echo "Nenhum arquivo informado. Encerrando."
  exit 1
fi

# Verificar se foram informados arquivos .html
for file in "${FILES_SELECTED[@]}"; do
  if [[ "$file" == *.html ]]; then
    echo "Erro: Arquivo '$file' é .html. MultiQC precisa dos arquivos .zip, não .html."
    echo "Encerrando."
    exit 1
  fi
  if [[ "$file" != *.zip ]]; then
    echo "Erro: Arquivo '$file' não é .zip. Apenas arquivos .zip são aceitos."
    echo "Encerrando."
    exit 1
  fi
  if [[ ! -f "$INPUT_DIR/$file" ]]; then
    echo "Erro: Arquivo '$file' não encontrado no diretório informado."
    echo "Encerrando."
    exit 1
  fi

  # Atualiza para usar caminho completo dentro do container
  FILES_TO_ANALYZE+=("/data/$file")

done

# Criar pasta de saída
TIMESTAMP=$(date +"%d-%m-%Y_%Hh%Mm")
OUTPUT_DIR="$INPUT_DIR/multiqc_output_$TIMESTAMP"
mkdir -p "$OUTPUT_DIR"

# Resumo para o usuário
echo "========= RESUMO FINAL ========="
echo "Diretório de entrada: $INPUT_DIR"
echo "Arquivos selecionados para MultiQC:"
for file in "${FILES_SELECTED[@]}"; do
  echo "- $file"
done
echo ""
echo "Diretório de saída: $OUTPUT_DIR"
echo "Docker UID: ${CUSTOM_UID}"
echo "================================="
echo ""
echo "Está tudo correto? (yes/no)"
read FINAL_CONFIRM
if [[ "$FINAL_CONFIRM" == "exit" ]]; then exit 0; fi
if [[ "$FINAL_CONFIRM" != "yes" ]]; then
  echo "Execução cancelada."
  exit 0
fi

# Rodar MultiQC
echo "Iniciando execução do MultiQC..."
docker run --rm \
  -v "$INPUT_DIR":/data \
  -v "$OUTPUT_DIR":/output \
  $USE_USER \
  ewels/multiqc \
  multiqc /data --outdir /output --force

# Finalização
echo ""
echo "==================================================="
echo "MultiQC finalizado."
echo "Relatório disponível em: $OUTPUT_DIR"
echo "==================================================="
exit 0
