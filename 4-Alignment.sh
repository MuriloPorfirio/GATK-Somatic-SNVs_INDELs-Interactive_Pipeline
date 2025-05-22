#!/bin/bash

# ===================== Introduction =====================
echo "==================================================="
echo " Paired-end FASTQ Alignment using BWA-MEM + SAM/BAM Conversion "
echo "==================================================="
echo ""
echo "This script aligns paired-end .fastq.gz files using the BWA-MEM algorithm"
echo "and converts the resulting SAM files into sorted BAM files."
echo "Required inputs:"
echo " - Paired FASTQ files (R1 and R2), concatenated and trimmed"
echo " - Reference genome index files (.fa, .amb, .ann, .bwt, .pac, .sa, .fai)"
echo "Outputs generated:"
echo " - Aligned and sorted BAM files"
echo " - BAM index files (.bai)"
echo " - Process log file"
echo ""
echo "IMPORTANT:"
echo " - Make sure your FASTQ files are in the correct format"
echo " - Use 'pwd' to copy paths correctly"
echo ""
echo "Can we continue? (yes/no)"
read START_CONFIRM
if [[ "$START_CONFIRM" != "yes" ]]; then
  echo "Script terminated."
  exit 0
fi

# ===================== Sample Inputs =====================
echo "Enter the FULL path to the directory containing the FASTQ files:"
echo "(Use 'pwd' to copy the path correctly, right-click to paste.)"
read INPUT_DIR
if [[ "$INPUT_DIR" == "exit" ]]; then exit 0; fi
if [[ ! -d "$INPUT_DIR" ]]; then
  echo "Error: Directory '$INPUT_DIR' not found."
  exit 1
fi

# List files
echo "Files found in the directory:"
ls "$INPUT_DIR"
echo ""
echo "Now enter the FASTQ file names you wish to align, separated by space."
echo "Example: patient1_R1.fastq.gz patient1_R2.fastq.gz patient2_R1.fastq.gz patient2_R2.fastq.gz"
read -a FASTQ_FILES
if [[ "$FASTQ_FILES" == "exit" ]]; then exit 0; fi

# Validate number of files
NUM_FILES=${#FASTQ_FILES[@]}
if (( $NUM_FILES % 2 != 0 )); then
  echo "Error: Odd number of files provided. It must be even (R1 and R2 for each sample)."
  exit 1
fi

# ===================== Reference Genome =====================
echo "Now enter the FULL path to the directory containing the reference genome and index files:"
read GENOME_DIR
if [[ "$GENOME_DIR" == "exit" ]]; then exit 0; fi
if [[ ! -d "$GENOME_DIR" ]]; then
  echo "Error: Directory '$GENOME_DIR' not found."
  exit 1
fi

echo "Files found in the genome directory:"
ls "$GENOME_DIR"
echo "Now enter the EXACT filenames for each required genome file:"
echo "(.fa, .amb, .ann, .bwt, .pac, .sa, .fai)"

read -p ".fa file: " GENOME_FA
read -p ".fa.amb file: " GENOME_AMB
read -p ".fa.ann file: " GENOME_ANN
read -p ".fa.bwt file: " GENOME_BWT
read -p ".fa.pac file: " GENOME_PAC
read -p ".fa.sa file: " GENOME_SA
read -p ".fa.fai file: " GENOME_FAI

# ===================== Execution Settings =====================
echo "Do you want to run Docker with a specific user (--user UID)? Enter the UID, or 'no' to run as root, or 'exit' to quit:"
read CUSTOM_UID
if [[ "$CUSTOM_UID" == "exit" ]]; then exit 0; fi
if [[ "$CUSTOM_UID" =~ ^[0-9]+$ ]]; then
  USE_USER="--user $CUSTOM_UID"
else
  USE_USER=""
fi

# Detect system resources
TOTAL_THREADS=$(nproc)
TOTAL_RAM=$(free -g | awk '/Mem:/ {print $2}')

echo "Your system has:"
echo "- Available threads: $TOTAL_THREADS"
echo "- Available RAM: ${TOTAL_RAM}GB"

echo "How many threads do you want to use for BWA? (recommended: leave a few free, e.g., $((TOTAL_THREADS - 4)))"
read NUM_THREADS

echo "How much RAM do you want to allocate? (min: 4GB, max: ${TOTAL_RAM}GB, recommend leaving a buffer)"
read MEMORY_RAM

# Create output directory
TIMESTAMP=$(date +"%d-%m-%Y_%Hh%Mm")
OUTPUT_DIR="$INPUT_DIR/4-aligned_BAM_files_$TIMESTAMP"
mkdir -p "$OUTPUT_DIR"

# ===================== Recap =====================
echo "========= RECAP ========="
echo "Input FASTQ files:"
for file in "${FASTQ_FILES[@]}"; do echo "- $file"; done
echo "FASTQ directory: $INPUT_DIR"
echo ""
echo "Reference genome files:"
echo "- $GENOME_FA"
echo "- $GENOME_AMB"
echo "- $GENOME_ANN"
echo "- $GENOME_BWT"
echo "- $GENOME_PAC"
echo "- $GENOME_SA"
echo "- $GENOME_FAI"
echo "Genome directory: $GENOME_DIR"
echo ""
echo "Threads to use: $NUM_THREADS"
echo "RAM to use: ${MEMORY_RAM}GB"
echo ""
echo "Outputs will be saved in: $OUTPUT_DIR"
echo "=========================="
echo "Shall we start? (yes/no)"
read FINAL_CONFIRM
if [[ "$FINAL_CONFIRM" != "yes" ]]; then
  echo "Execution cancelled."
  exit 0
fi

# ===================== Execution =====================
echo "Starting alignment and conversion..."

SCRIPT_START_TIME=$(date +%s)

for ((i=0; i<${#FASTQ_FILES[@]}; i+=2)); do
  SAMPLE_R1=${FASTQ_FILES[$i]}
  SAMPLE_R2=${FASTQ_FILES[$((i+1))]}
  SAMPLE_NAME=$(basename "$SAMPLE_R1" | cut -d '_' -f1)
  OUTPUT_SAM="$OUTPUT_DIR/${SAMPLE_NAME}.sam"
  OUTPUT_BAM="$OUTPUT_DIR/${SAMPLE_NAME}_aligned.bam"
  LOG_FILE="$OUTPUT_DIR/${SAMPLE_NAME}_alignment.log"

  SAMPLE_START_TIME=$(date +%s)

  # Run BWA to generate SAM
  docker run --rm \
    $USE_USER \
    -v "$INPUT_DIR":/data \
    -v "$GENOME_DIR":/ref \
    -v "$OUTPUT_DIR":/output \
    biocontainers/bwa:v0.7.17_cv1 \
    bash -c "set -e; bwa mem -t $NUM_THREADS /ref/$GENOME_FA /data/$SAMPLE_R1 /data/$SAMPLE_R2 > /output/${SAMPLE_NAME}.sam" 2>&1 | tee "$LOG_FILE"

  # Run Samtools to convert SAM to sorted BAM
  docker run --rm \
    $USE_USER \
    -v "$OUTPUT_DIR":/workdir \
    staphb/samtools:latest \
    bash -c "set -e; samtools view -@ $NUM_THREADS -bS /workdir/${SAMPLE_NAME}.sam | samtools sort -@ $NUM_THREADS -o /workdir/${SAMPLE_NAME}_aligned.bam -; samtools index /workdir/${SAMPLE_NAME}_aligned.bam"

  # Remove intermediate SAM file
  rm "$OUTPUT_DIR/${SAMPLE_NAME}.sam"

  SAMPLE_END_TIME=$(date +%s)
  SAMPLE_DURATION=$((SAMPLE_END_TIME - SAMPLE_START_TIME))
  echo "Time to process ${SAMPLE_NAME}: ${SAMPLE_DURATION} seconds." | tee -a "$LOG_FILE"
done

SCRIPT_END_TIME=$(date +%s)
TOTAL_DURATION=$((SCRIPT_END_TIME - SCRIPT_START_TIME))
TOTAL_HOURS=$((TOTAL_DURATION / 3600))
TOTAL_MINUTES=$(((TOTAL_DURATION % 3600) / 60))
TOTAL_SECONDS=$((TOTAL_DURATION % 60))

echo ""
echo "======================================================"
echo " Alignment completed! BAM files and logs saved in:"
echo " $OUTPUT_DIR"
echo " Total execution time: ${TOTAL_HOURS}h ${TOTAL_MINUTES}m ${TOTAL_SECONDS}s"
echo "======================================================"
exit 0
