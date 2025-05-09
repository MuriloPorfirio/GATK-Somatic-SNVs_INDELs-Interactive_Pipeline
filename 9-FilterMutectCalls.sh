#!/bin/bash

echo "=============================================================="
echo "SCRIPT INTERATIVO PARA FILTRAGEM DE VARIANTES SOMÁTICAS"
echo "ETAPA: FilterMutectCalls (GATK)"
echo "=============================================================="
echo ""
echo "Esta etapa aplica filtros estatísticos e heurísticos às variantes somáticas"
echo "detectadas pelo Mutect2. O objetivo é eliminar falsos positivos e produzir"
echo "um conjunto confiável de variantes para downstream (ex: anotações, validações, etc)."
echo ""
echo "FERRAMENTA UTILIZADA:"
echo "- GATK FilterMutectCalls (via Docker: broadinstitute/gatk)"
echo ""
echo "INPUTS NECESSÁRIOS:"
echo "- Arquivo .vcf.gz gerado pelo Mutect2"
echo "- Genoma de referência (.fa) com .fai e .dict na mesma pasta"
echo ""
echo "SAÍDA:"
echo "- VCF filtrado (.vcf.gz)"
echo "- Índice (.vcf.gz.tbi)"
echo "- Arquivo de estatísticas"
echo ""

echo "Podemos continuar? (yes/no)"
read START_CONFIRM
[[ "$START_CONFIRM" != "yes" ]] && echo "Execução cancelada." && exit 0

# === Genoma ===
echo ""
echo "Digite o caminho COMPLETO da pasta com o genoma .fa:"
read GENOME_DIR
[[ ! -d "$GENOME_DIR" ]] && echo "Erro: pasta não encontrada." && exit 1
echo "Arquivos encontrados:"
ls "$GENOME_DIR"
echo "Digite o NOME do arquivo .fa de referência (ex: hg38.fa):"
read GENOME_FA
GENOME="$GENOME_DIR/$GENOME_FA"

if [[ ! -f "$GENOME" || ! -f "$GENOME_DIR/${GENOME_FA}.fai" || ! -f "$GENOME_DIR/${GENOME_FA%.fa}.dict" ]]; then
  echo "Erro: .fa, .fai e .dict devem estar todos na mesma pasta com nomes consistentes."
  exit 1
fi

# === VCF do Mutect2 ===
echo ""
echo "Digite o caminho COMPLETO da pasta com o VCF gerado pelo Mutect2:"
read VCF_DIR
[[ ! -d "$VCF_DIR" ]] && echo "Erro: pasta não encontrada." && exit 1
echo "Arquivos encontrados:"
ls "$VCF_DIR"
echo "Digite o NOME EXATO do arquivo .vcf.gz do Mutect2:"
read VCF_MUTECT
[[ ! -f "$VCF_DIR/$VCF_MUTECT" ]] && echo "Erro: VCF não encontrado." && exit 1
[[ ! -f "$VCF_DIR/${VCF_MUTECT}.tbi" ]] && echo "Erro: arquivo .tbi correspondente não encontrado." && exit 1

# === UID Docker ===
echo ""
echo "Digite o UID para rodar no Docker (ex: 1006) ou 'no' para root:"
read UID_INPUT
USE_USER=""
[[ "$UID_INPUT" =~ ^[0-9]+$ ]] && USE_USER="--user $UID_INPUT"

# === Threads e RAM ===
TOTAL_THREADS=$(nproc)
echo "Seu sistema possui $TOTAL_THREADS threads. Quantas deseja usar?"
read THREADS
TOTAL_RAM=$(free -g | awk '/Mem:/ {print $2}')
echo "Memória total: ${TOTAL_RAM} GB. Quanto deseja usar?"
read RAM

# === Criar pasta de saída ===
TS=$(date +"%d-%m-%Y_%Hh%Mm")
OUTPUT_DIR="$VCF_DIR/9-FilterMutectCalls_output_$TS"
mkdir -p "$OUTPUT_DIR"

# === Recapitulação ===
echo ""
echo "================ RECAPITULAÇÃO ================"
echo "Genoma: $GENOME"
echo "VCF Mutect2: $VCF_MUTECT"
echo "Pasta de saída: $OUTPUT_DIR"
echo "Usuário Docker: ${UID_INPUT:-root}"
echo "Threads: $THREADS | RAM: ${RAM}GB"
echo "================================================"
echo "Confirmar execução? (yes/no)"
read CONFIRM
[[ "$CONFIRM" != "yes" ]] && echo "Execução cancelada." && exit 0

# === Execução ===
VCF_BASENAME=$(basename "$VCF_MUTECT" .vcf.gz)
VCF_IN="/vcf/$VCF_MUTECT"
VCF_OUT="${VCF_BASENAME}_filtered.vcf.gz"
LOG_OUT="${VCF_BASENAME}_filtermutectcalls.log"

SCRIPT_START=$(date +%s)

docker run --rm $USE_USER \
  -v "$GENOME_DIR":/genome \
  -v "$VCF_DIR":/vcf \
  -v "$OUTPUT_DIR":/output \
  broadinstitute/gatk \
  bash -c "gatk FilterMutectCalls \
    -R /genome/$GENOME_FA \
    -V $VCF_IN \
    -O /output/$VCF_OUT \
    --verbosity INFO" \
  &> "$OUTPUT_DIR/$LOG_OUT"

SCRIPT_END=$(date +%s)
DURATION=$((SCRIPT_END - SCRIPT_START))

echo ""
echo "=============================================================="
echo "FilterMutectCalls finalizado com sucesso!"
echo "Saída disponível em: $OUTPUT_DIR"
echo "Tempo total: $((DURATION / 3600))h $((DURATION % 3600 / 60))m $((DURATION % 60))s"
echo "=============================================================="
