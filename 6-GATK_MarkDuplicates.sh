#!/bin/bash

# ===================== Apresentação =====================
echo "=============================================================="
echo "        Script interativo para GATK MarkDuplicates"
echo "=============================================================="
echo ""
echo "O GATK MarkDuplicates identifica e marca leituras duplicadas geradas por PCR"
echo "para evitar vieses na análise de variantes."
echo ""
echo "IMPORTÂNCIA DA FERRAMENTA:"
echo " - Remove vieses causados por duplicatas PCR."
echo " - Melhora a precisão das variantes chamadas."
echo " - Mantém profundidade real (sem profundidade artificial)."
echo ""
echo "REQUISITOS ANTES DE EXECUTAR ESSE SCRIPT:"
echo " - Arquivos BAM já alinhados e ordenados (gerados após alinhamento com BWA e conversão com Samtools)."
echo ""
echo "Podemos continuar? (yes/no)"
read START_CONFIRM
if [[ "$START_CONFIRM" != "yes" ]]; then
  echo "Script encerrado."
  exit 0
fi

# ===================== Input dos arquivos BAM =====================
echo "Digite o caminho COMPLETO onde estão os arquivos BAM alinhados:"
echo "(Use 'pwd' para copiar corretamente, cole com botão direito.)"
read INPUT_DIR
if [[ "$INPUT_DIR" == "exit" ]]; then exit 0; fi
if [[ ! -d "$INPUT_DIR" ]]; then
  echo "Erro: Diretório '$INPUT_DIR' não encontrado."
  exit 1
fi

echo "Arquivos encontrados no diretório escolhido:"
ls "$INPUT_DIR"
echo ""
echo "Digite exatamente o nome dos arquivos BAM que você deseja processar, separados por espaço."
echo "Exemplo: paciente1.bam paciente2.bam paciente3.bam"
read -a BAM_FILES
if [[ "$BAM_FILES" == "exit" ]]; then exit 0; fi

# ===================== Usuário Docker =====================
echo "Deseja rodar o Docker com um usuário específico (--user UID)?"
echo "Informe o UID (por exemplo: 1006), ou 'no' para rodar como root, ou 'exit' para sair:"
read CUSTOM_UID
if [[ "$CUSTOM_UID" == "exit" ]]; then exit 0; fi
if [[ "$CUSTOM_UID" =~ ^[0-9]+$ ]]; then
  USE_USER="--user $CUSTOM_UID"
else
  USE_USER=""
fi

# ===================== Recursos Computacionais =====================
echo "O GATK MarkDuplicates usa CPU. O padrão do GATK é usar 1 thread."
TOTAL_THREADS=$(nproc)
echo "Seu servidor possui $TOTAL_THREADS threads disponíveis."
echo "Quantas threads você quer utilizar? (recomendado: deixar algumas livres, ex: $((TOTAL_THREADS - 2)))"
read NUM_THREADS

DEFAULT_RAM=8
echo "O GATK MarkDuplicates por padrão usa em torno de ${DEFAULT_RAM}GB de RAM."
TOTAL_RAM=$(free -g | awk '/Mem:/ {print $2}')
echo "Seu servidor tem ${TOTAL_RAM}GB de RAM."
echo "Quantos GB de RAM você deseja utilizar? (exemplo recomendado: $((TOTAL_RAM - 2)))"
read MEMORY_RAM

# ===================== Criando diretório de output =====================
TIMESTAMP=$(date +"%d-%m-%Y_%Hh%Mm")
OUTPUT_DIR="$INPUT_DIR/6-output_marked_duplicates_$TIMESTAMP"
mkdir -p "$OUTPUT_DIR"

# ===================== Recapitulação =====================
echo "================ RECAPITULAÇÃO ================"
echo "Arquivos selecionados:"
for file in "${BAM_FILES[@]}"; do echo "- $file"; done
echo ""
echo "Diretório dos arquivos BAM: $INPUT_DIR"
echo "Usuário Docker: ${CUSTOM_UID:-root}"
echo "Threads escolhidas: $NUM_THREADS"
echo "Memória RAM escolhida: ${MEMORY_RAM}GB"
echo "Outputs serão salvos em: $OUTPUT_DIR"
echo "================================================="
echo "Podemos começar? (yes/no)"
read FINAL_CONFIRM
if [[ "$FINAL_CONFIRM" != "yes" ]]; then
  echo "Execução cancelada."
  exit 0
fi

# ===================== Execução =====================
echo "Iniciando a execução do GATK MarkDuplicates..."
SCRIPT_START_TIME=$(date +%s)

for BAM_FILE in "${BAM_FILES[@]}"; do
  SAMPLE_NAME=$(basename "$BAM_FILE" .bam)
  OUTPUT_BAM="$OUTPUT_DIR/${SAMPLE_NAME}_marked_duplicates.bam"
  METRICS_FILE="$OUTPUT_DIR/${SAMPLE_NAME}_duplication_metrics.txt"
  LOG_FILE="$OUTPUT_DIR/${SAMPLE_NAME}_markduplicates.log"

  SAMPLE_START_TIME=$(date +%s)

  docker run --rm $USE_USER \
    -v "$INPUT_DIR":/data \
    -v "$OUTPUT_DIR":/output \
    broadinstitute/gatk:latest gatk MarkDuplicates \
    -I /data/$BAM_FILE \
    -O /output/$(basename $OUTPUT_BAM) \
    -M /output/$(basename $METRICS_FILE) \
    --QUIET false \
    --VERBOSITY INFO \
    --TMP_DIR /output/tmp_markdup \
    --MAX_FILE_HANDLES_FOR_READ_ENDS_MAP 8000 \
    --REMOVE_DUPLICATES false \
    --ASSUME_SORTED true \
    --CREATE_INDEX true \
    --VALIDATION_STRINGENCY LENIENT \
    --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 \
    2>&1 | tee "$LOG_FILE"

  SAMPLE_END_TIME=$(date +%s)
  SAMPLE_DURATION=$((SAMPLE_END_TIME - SAMPLE_START_TIME))
  echo "Tempo para processar ${SAMPLE_NAME}: ${SAMPLE_DURATION} segundos." | tee -a "$LOG_FILE"
done

SCRIPT_END_TIME=$(date +%s)
TOTAL_DURATION=$((SCRIPT_END_TIME - SCRIPT_START_TIME))
TOTAL_HOURS=$((TOTAL_DURATION / 3600))
TOTAL_MINUTES=$(((TOTAL_DURATION % 3600) / 60))
TOTAL_SECONDS=$((TOTAL_DURATION % 60))

echo ""
echo "=============================================================="
echo "GATK MarkDuplicates concluído! Arquivos BAM e logs estão em:"
echo "$OUTPUT_DIR"
echo "Tempo total de execução: ${TOTAL_HOURS}h ${TOTAL_MINUTES}m ${TOTAL_SECONDS}s."
echo "=============================================================="
exit 0
