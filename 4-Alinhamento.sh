#!/bin/bash

# ===================== Apresentacao =====================
echo "==================================================="
echo " Script de Alinhamento de FASTQ pareados usando BWA-MEM "
echo "==================================================="
echo ""
echo "Este script alinha arquivos FASTQ.gz pareados usando o algoritmo BWA-MEM."
echo "Inputs necessários:"
echo " - Arquivos FASTQ pareados (R1 e R2), concatenados e trimados."
echo " - Arquivos de índice do genoma de referência (.fa, .amb, .ann, .bwt, .pac, .sa, .fai)."
echo "Outputs gerados:"
echo " - Arquivos BAM alinhados e ordenados."
echo " - Um arquivo de LOG do processo."
echo ""
echo "IMPORTANTE:"
echo " - Tenha certeza de que seus arquivos FASTQ estão no formato correto."
echo " - Use o comando 'pwd' para copiar caminhos corretamente."
echo ""
echo "Podemos continuar? (yes/no)"
read START_CONFIRM
if [[ "$START_CONFIRM" != "yes" ]]; then
  echo "Script encerrado."
  exit 0
fi

# ===================== Inputs de amostras =====================
echo "Digite o caminho COMPLETO onde estão os arquivos FASTQ:"
echo "(Use 'pwd' para copiar corretamente, cole com botão direito.)"
read INPUT_DIR
if [[ "$INPUT_DIR" == "exit" ]]; then exit 0; fi
if [[ ! -d "$INPUT_DIR" ]]; then
  echo "Erro: Diretório '$INPUT_DIR' não encontrado."
  exit 1
fi

# Listar arquivos
echo "Arquivos encontrados no diretório:"
ls "$INPUT_DIR"
echo ""
echo "Agora informe os arquivos FASTQ que deseja alinhar, separados por espaço."
echo "Exemplo: paciente1_R1.fastq.gz paciente1_R2.fastq.gz paciente2_R1.fastq.gz paciente2_R2.fastq.gz"
read -a FASTQ_FILES
if [[ "$FASTQ_FILES" == "exit" ]]; then exit 0; fi

# Validar número de arquivos (pares)
NUM_FILES=${#FASTQ_FILES[@]}
if (( $NUM_FILES % 2 != 0 )); then
  echo "Erro: número ímpar de arquivos informados. Deve ser par (R1 e R2 para cada amostra)."
  exit 1
fi

# ===================== Inputs de genoma =====================
echo "Agora digite o caminho COMPLETO do diretório onde está o genoma de referência e seus índices:"
read GENOME_DIR
if [[ "$GENOME_DIR" == "exit" ]]; then exit 0; fi
if [[ ! -d "$GENOME_DIR" ]]; then
  echo "Erro: Diretório '$GENOME_DIR' não encontrado."
  exit 1
fi

echo "Arquivos encontrados no diretório de genoma:"
ls "$GENOME_DIR"
echo "Agora informe exatamente os nomes dos arquivos de genoma para cada tipo exigido:"
echo "(.fa, .amb, .ann, .bwt, .pac, .sa, .fai)"

read -p ".fa arquivo: " GENOME_FA
read -p ".fa.amb arquivo: " GENOME_AMB
read -p ".fa.ann arquivo: " GENOME_ANN
read -p ".fa.bwt arquivo: " GENOME_BWT
read -p ".fa.pac arquivo: " GENOME_PAC
read -p ".fa.sa arquivo: " GENOME_SA
read -p ".fa.fai arquivo: " GENOME_FAI

# ===================== Configurações de Execução =====================
echo "Deseja rodar o Docker com um usuário específico (--user UID)? Informe o UID, ou 'no' para rodar como root, ou 'exit' para sair:"
read CUSTOM_UID
if [[ "$CUSTOM_UID" == "exit" ]]; then exit 0; fi
if [[ "$CUSTOM_UID" =~ ^[0-9]+$ ]]; then
  USE_USER="--user $CUSTOM_UID"
else
  USE_USER=""
fi

# Descobrir número de threads e RAM
TOTAL_THREADS=$(nproc)
TOTAL_RAM=$(free -g | awk '/Mem:/ {print $2}')

echo "Seu servidor possui:"
echo "- Threads disponíveis: $TOTAL_THREADS"
echo "- Memória RAM disponível: ${TOTAL_RAM}GB"

echo "Quantas threads deseja usar? (recomendo deixar um pouco livre, usar tipo $((TOTAL_THREADS-4)))"
read NUM_THREADS

echo "Quanto de RAM deseja alocar? (min: 4GB, máx: ${TOTAL_RAM}GB, recomendo deixar um pouco livre)"
read MEMORY_RAM

# Criar pasta de output
TIMESTAMP=$(date +"%d-%m-%Y_%Hh%Mm")
OUTPUT_DIR="$INPUT_DIR/arquivos_alinhados_(BAM)_$TIMESTAMP"
mkdir -p "$OUTPUT_DIR"

# ===================== Recapitulação =====================
echo "========= RECAPITULAÇÃO ========="
echo "Inputs das amostras:"
for file in "${FASTQ_FILES[@]}"; do echo "- $file"; done
echo "Diretório das amostras: $INPUT_DIR"
echo ""
echo "Inputs do genoma de referência:"
echo "- $GENOME_FA"
echo "- $GENOME_AMB"
echo "- $GENOME_ANN"
echo "- $GENOME_BWT"
echo "- $GENOME_PAC"
echo "- $GENOME_SA"
echo "- $GENOME_FAI"
echo "Diretório do genoma: $GENOME_DIR"
echo ""
echo "Threads a utilizar: $NUM_THREADS"
echo "RAM a utilizar: ${MEMORY_RAM}GB"
echo ""
echo "Outputs serão salvos em: $OUTPUT_DIR"
echo "====================================="
echo "Podemos começar? (yes/no)"
read FINAL_CONFIRM
if [[ "$FINAL_CONFIRM" != "yes" ]]; then
  echo "Execução cancelada."
  exit 0
fi

# ===================== Execução =====================
echo "Iniciando o alinhamento..."

for ((i=0; i<${#FASTQ_FILES[@]}; i+=2)); do
  SAMPLE_R1=${FASTQ_FILES[$i]}
  SAMPLE_R2=${FASTQ_FILES[$((i+1))]}
  SAMPLE_NAME=$(basename "$SAMPLE_R1" | cut -d '_' -f1)
  OUTPUT_BAM="$OUTPUT_DIR/${SAMPLE_NAME}_aligned.bam"
  LOG_FILE="$OUTPUT_DIR/${SAMPLE_NAME}_alignment.log"

  docker run --rm \
    $USE_USER \
    -v "$INPUT_DIR":/data \
    -v "$GENOME_DIR":/ref \
    -v "$OUTPUT_DIR":/output \
    biocontainers/bwa:v0.7.17_cv1 \
    bash -c "set -e; bwa mem -t $NUM_THREADS /ref/$GENOME_FA /data/$SAMPLE_R1 /data/$SAMPLE_R2 | samtools view -Sb - | samtools sort -o /output/${SAMPLE_NAME}_aligned.bam -" 2>&1 | tee "$LOG_FILE"

done

echo ""
echo "======================================================"
echo " Alinhamento concluído! Arquivos BAM e logs salvos em:"
echo " $OUTPUT_DIR"
echo "======================================================"
exit 0
