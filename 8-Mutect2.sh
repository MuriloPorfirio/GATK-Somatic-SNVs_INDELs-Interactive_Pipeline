#!/bin/bash

echo "=============================================================="
echo " INTERACTIVE SCRIPT FOR SOMATIC VARIANT CALLING (Mutect2)"
echo "=============================================================="
echo ""
echo "This script runs GATK Mutect2 to identify somatic variants"
echo "in cancer samples. It supports both TUMOR-ONLY and TUMOR/NORMAL modes."
echo ""
echo "TOOLS USED:"
echo "- GATK (via Docker: broadinstitute/gatk)"
echo ""
echo "THIS SCRIPT WILL GENERATE:"
echo "- A .vcf.gz file with somatic variants"
echo "- A .vcf.gz.tbi index file"
echo "- A log file of the execution"
echo ""
echo "REQUIRED FILES:"
echo "- Tumor BAM file (with RG and marked duplicates)"
echo "- Reference genome .fa (.fai and .dict in the same folder)"
echo ""
echo "OPTIONAL (highly recommended):"
echo "- Panel of Normals (PoN): recurrent artifacts in normal tissues"
echo "- Matched Normal BAM (germline tissue from the same patient)"
echo "- Germline Resource: gnomAD allele frequency file (e.g., af-only-gnomad.hg38.vcf.gz + .tbi)"
echo ""
echo "NOTE: If you don’t have a PoN or normal sample, using gnomAD is highly recommended."
echo "You can download gnomAD from: https://console.cloud.google.com/storage/browser/gatk-best-practices/somatic-hg38"
echo ""

echo "Shall we continue? (yes/no)"
read START_CONFIRM
[[ "$START_CONFIRM" != "yes" ]] && echo "Script canceled." && exit 0

# === Reference Genome ===
echo ""
echo "Enter the FULL path to the folder containing the .fa genome file:"
read GENOME_DIR
[[ ! -d "$GENOME_DIR" ]] && echo "Error: directory not found." && exit 1
echo "Files found:"
ls "$GENOME_DIR"
echo "Enter the EXACT NAME of the .fa file:"
read GENOME_FA
GENOME="$GENOME_DIR/$GENOME_FA"

if [[ ! -f "$GENOME" || ! -f "$GENOME_DIR/${GENOME_FA}.fai" || ! -f "$GENOME_DIR/${GENOME_FA%.fa}.dict" ]]; then
  echo "Error: .fa, .fai, and .dict must all be present and consistent."
  exit 1
fi

# === Tumor BAM ===
echo ""
echo "Enter the FULL path to the folder with the BAM files:"
read BAM_DIR
[[ ! -d "$BAM_DIR" ]] && echo "Error: directory not found." && exit 1
echo "Files found:"
ls "$BAM_DIR"
echo "Enter the EXACT NAME of the TUMOR BAM file:"
read BAM_TUMOR
[[ ! -f "$BAM_DIR/$BAM_TUMOR" ]] && echo "Error: file not found." && exit 1

# === Panel of Normals ===
echo ""
echo "Do you have a Panel of Normals (PoN)? (yes/no)"
read HAS_PON
if [[ "$HAS_PON" == "yes" ]]; then
  echo "Enter the FULL path to the folder with the PoN .vcf.gz file:"
  read PON_DIR
  [[ ! -d "$PON_DIR" ]] && echo "Error: directory not found." && exit 1
  echo "Files found:"
  ls "$PON_DIR"
  echo "Enter the EXACT NAME of the PoN .vcf.gz file:"
  read PON_FILE
  PON_PATH="$PON_DIR/$PON_FILE"
  [[ ! -f "$PON_PATH" ]] && echo "Error: PoN file not found." && exit 1
  if [[ ! -f "$PON_PATH.tbi" ]]; then
    echo "Error: corresponding .tbi index file for the PoN is missing."
    exit 1
  fi
fi

# === Matched Normal BAM ===
echo ""
echo "Do you have a matched NORMAL BAM file? (yes/no)"
read HAS_NORMAL
if [[ "$HAS_NORMAL" == "yes" ]]; then
  echo "Enter the FULL path to the folder with the NORMAL BAM file:"
  read NORMAL_DIR
  [[ ! -d "$NORMAL_DIR" ]] && echo "Error: directory not found." && exit 1
  echo "Files found:"
  ls "$NORMAL_DIR"
  echo "Enter the EXACT NAME of the NORMAL BAM file:"
  read BAM_NORMAL
  BAM_NORMAL_PATH="$NORMAL_DIR/$BAM_NORMAL"
  [[ ! -f "$BAM_NORMAL_PATH" ]] && echo "Error: normal BAM file not found." && exit 1
  if [[ ! -f "$BAM_NORMAL_PATH.bai" ]]; then
    echo "Error: .bai index file for the normal BAM is missing."
    exit 1
  fi
fi

# === Germline Resource ===
echo ""
echo "Do you have the af-only-gnomad.hg38.vcf.gz (gnomAD)? (yes/no)"
read HAS_GNOMAD
if [[ "$HAS_GNOMAD" == "yes" ]]; then
  echo "Enter the FULL path to the folder containing the gnomAD file:"
  read GNOMAD_DIR
  [[ ! -d "$GNOMAD_DIR" ]] && echo "Error: directory not found." && exit 1
  echo "Files found:"
  ls "$GNOMAD_DIR"
  echo "Enter the EXACT NAME of the gnomAD VCF file:"
  read GNOMAD_VCF
  GNOMAD_VCF_PATH="$GNOMAD_DIR/$GNOMAD_VCF"
  [[ ! -f "$GNOMAD_VCF_PATH" ]] && echo "Error: gnomAD VCF not found." && exit 1
  if [[ ! -f "$GNOMAD_VCF_PATH.tbi" ]]; then
    echo "Error: .tbi index for gnomAD VCF not found."
    exit 1
  fi
else
  echo "WARNING: Without gnomAD, Mutect2 may have reduced ability to filter germline variants."
fi

# === Docker UID ===
echo ""
echo "Enter the UID to run Docker (e.g., 1006), or 'no' to run as root:"
read UID_INPUT
USE_USER=""
[[ "$UID_INPUT" =~ ^[0-9]+$ ]] && USE_USER="--user $UID_INPUT"

# === Threads and RAM ===
TOTAL_THREADS=$(nproc)
echo "Your system has $TOTAL_THREADS available threads. How many do you want to use?"
read THREADS
TOTAL_RAM=$(free -g | awk '/Mem:/ {print $2}')
echo "Total memory: ${TOTAL_RAM} GB. How much do you want to allocate?"
read RAM

# === Output Directory ===
TS=$(date +"%d-%m-%Y_%Hh%Mm")
OUTPUT_DIR="$BAM_DIR/8-mutect2_output_$TS"
mkdir -p "$OUTPUT_DIR"

# === Summary ===
echo ""
echo "================ SUMMARY ================"
echo "Reference genome: $GENOME"
echo "Tumor BAM: $BAM_TUMOR"
[[ "$HAS_NORMAL" == "yes" ]] && echo "Normal BAM: $BAM_NORMAL"
[[ "$HAS_PON" == "yes" ]] && echo "PoN: $PON_FILE"
[[ "$HAS_GNOMAD" == "yes" ]] && echo "gnomAD: $GNOMAD_VCF"
echo "Output directory: $OUTPUT_DIR"
echo "Docker user: ${UID_INPUT:-root}"
echo "Threads: $THREADS | RAM: ${RAM}GB"
echo "========================================="
echo "Confirm execution? (yes/no)"
read CONFIRM
[[ "$CONFIRM" != "yes" ]] && echo "Execution canceled." && exit 0

# === Execution ===
BAM_TUMOR_BASENAME=$(basename "$BAM_TUMOR" .bam)
VCF_OUT="$OUTPUT_DIR/${BAM_TUMOR_BASENAME}_somatic.vcf.gz"
LOG="$OUTPUT_DIR/${BAM_TUMOR_BASENAME}_mutect2.log"

CMD="gatk Mutect2 -R /genome/$(basename "$GENOME") -I /bam/$(basename "$BAM_TUMOR") -O /output/$(basename "$VCF_OUT")"

[[ "$HAS_NORMAL" == "yes" ]] && CMD+=" -I /normal/$(basename "$BAM_NORMAL") --normal-sample NORMAL"
[[ "$HAS_PON" == "yes" ]] && CMD+=" --panel-of-normals /pon/$(basename "$PON_FILE")"
[[ "$HAS_GNOMAD" == "yes" ]] && CMD+=" --germline-resource /gnomad/$(basename "$GNOMAD_VCF")"

docker run --rm $USE_USER \
  -v "$GENOME_DIR":/genome \
  -v "$BAM_DIR":/bam \
  $( [[ "$HAS_NORMAL" == "yes" ]] && echo -v "$NORMAL_DIR":/normal ) \
  $( [[ "$HAS_PON" == "yes" ]] && echo -v "$PON_DIR":/pon ) \
  $( [[ "$HAS_GNOMAD" == "yes" ]] && echo -v "$GNOMAD_DIR":/gnomad ) \
  -v "$OUTPUT_DIR":/output \
  broadinstitute/gatk \
  bash -c "$CMD" &> "$LOG"

echo "=============================================================="
echo "Mutect2 completed. Output available at: $OUTPUT_DIR"
echo "Log: $LOG"
echo "=============================================================="
