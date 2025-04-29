#!/bin/bash

# ===================== Apresentação =====================
echo "=============================================================="
echo "    Script interativo para GATK AddOrReplaceReadGroups"
echo "=============================================================="
echo ""
echo "Esta ferramenta adiciona cabeçalhos @RG (Read Group) a arquivos BAM,"
echo "pré-requisito para uso com GATK MarkDuplicates e outras ferramentas."
echo ""
echo "IMPORTÂNCIA DA FERRAMENTA:"
echo " - Permite que o GATK identifique grupos de leitura distintos."
echo " - Evita erros como NullPointerException ao rodar MarkDuplicates."
echo ""
echo "REQUISITOS: BAM alinhado, ordenado e sem read groups."
echo ""
echo "Podemos continuar? (yes/no)"
read START_CONFIRM
if [[ "$START_CONFIRM" != "yes" ]]; then
  echo "Script encerrado."
  exit 0
fi

# ===================== Diretório e arquivos =====================
echo "Digite o caminho COMPLETO onde está o arquivo BAM sem read groups:"
read INPUT_DIR
if [[ ! -d "$INPUT_DIR" ]]; then
  echo "Erro: diretório não encontrado."
  exit 1
fi

echo "Arquivos encontrados no diretório:"
ls "$INPUT_DIR"
echo "Digite o nome exato do arquivo BAM a ser processado:"
read BAM_FILE

# ===================== Identificação do Read Group =====================
echo "Agora insira as informações necessárias para identificar corretamente a amostra no GATK."
echo "Esses valores serão usados no cabeçalho do arquivo BAM."
echo ""
echo "RGID: ID do grupo de leitura. Pode ser o nome do arquivo ou da amostra."
read -p "Digite o RGID (ex: paciente01): " RGID

echo "RGLB: Identificador da biblioteca de preparo. Se você não sabe, use algo simples como 'lib1'."
read -p "Digite o RGLB (ex: lib1): " RGLB

echo "RGPL: Plataforma de sequenciamento usada. Escolha uma das opções abaixo:"
echo "1 = ILLUMINA"
echo "2 = IONTORRENT"
echo "3 = PACBIO"
echo "4 = ONT"
echo "5 = BGI"
read -p "Escolha o número correspondente: " RGPL_CHOICE
case $RGPL_CHOICE in
  1) RGPL="ILLUMINA";;
  2) RGPL="IONTORRENT";;
  3) RGPL="PACBIO";;
  4) RGPL="ONT";;
  5) RGPL="BGI";;
  *) echo "Opção inválida."; exit 1;;
esac

echo "RGPU: Unidade da plataforma. Pode inventar algo simples, como 'unit1'."
read -p "Digite o RGPU: " RGPU

echo "RGSM: Nome da amostra. Use o mesmo nome que você usará em todo o pipeline."
read -p "Digite o RGSM (ex: paciente01): " RGSM

# ===================== Usuário Docker =====================
echo "Deseja rodar o Docker com um usuário específico (--user UID)?"
echo "Informe o UID (exemplo: 1006), ou 'no' para rodar como root."
read CUSTOM_UID
if [[ "$CUSTOM_UID" =~ ^[0-9]+$ ]]; then
  USE_USER="--user $CUSTOM_UID"
else
  USE_USER=""
fi

# ===================== Diretório de saída =====================
TIMESTAMP=$(date +"%d-%m-%Y_%Hh%Mm")
OUTPUT_DIR="$INPUT_DIR/4-output_add_readgroups_$TIMESTAMP"
mkdir -p "$OUTPUT_DIR"
OUTPUT_BAM="$OUTPUT_DIR/${BAM_FILE%.bam}_with_RG.bam"
LOG_FILE="$OUTPUT_DIR/${BAM_FILE%.bam}_addRG.log"

# ===================== Recapitulação =====================
echo "================ RECAPITULAÇÃO ================"
echo "Arquivo de entrada: $BAM_FILE"
echo "Diretório de entrada: $INPUT_DIR"
echo "RGID: $RGID"
echo "RGLB: $RGLB"
echo "RGPL: $RGPL"
echo "RGPU: $RGPU"
echo "RGSM: $RGSM"
echo "Usuário Docker: ${CUSTOM_UID:-root}"
echo "Saída: $OUTPUT_BAM"
echo "Log: $LOG_FILE"
echo "================================================="
echo "Podemos começar? (yes/no)"
read FINAL_CONFIRM
if [[ "$FINAL_CONFIRM" != "yes" ]]; then
  echo "Execução cancelada."
  exit 0
fi

# ===================== Execução =====================
echo "Iniciando inserção de Read Groups..."
START_TIME=$(date +%s)

docker run --rm $USE_USER \
  -v "$INPUT_DIR":/data \
  -v "$OUTPUT_DIR":/output \
  broadinstitute/gatk:latest gatk AddOrReplaceReadGroups \
  -I /data/$BAM_FILE \
  -O /output/$(basename $OUTPUT_BAM) \
  -RGID $RGID \
  -RGLB $RGLB \
  -RGPL $RGPL \
  -RGPU $RGPU \
  -RGSM $RGSM \
  2>&1 | tee "$LOG_FILE"

END_TIME=$(date +%s)
DURATION=$((END_TIME - START_TIME))
HOURS=$((DURATION / 3600))
MINUTES=$(((DURATION % 3600) / 60))
SECONDS=$((DURATION % 60))

echo ""
echo "=============================================================="
echo "Processo finalizado! BAM com read groups salvo em:"
echo "$OUTPUT_BAM"
echo "Tempo total: ${HOURS}h ${MINUTES}m ${SECONDS}s"
echo "=============================================================="
exit 0
