#!/bin/bash

# ===================== Presentation =====================
echo "=============================================================="
echo "    Interactive Script for GATK AddOrReplaceReadGroups"
echo "=============================================================="
echo ""
echo "This tool adds @RG (Read Group) headers to BAM files,"
echo "which is a prerequisite for using GATK MarkDuplicates and other tools."
echo ""
echo "IMPORTANCE OF THIS TOOL:"
echo " - Allows GATK to distinguish between different read groups."
echo " - Prevents errors such as NullPointerException when running MarkDuplicates."
echo ""
echo "REQUIREMENTS: Aligned and sorted BAM without existing read groups."
echo ""
echo "Shall we continue? (yes/no)"
read START_CONFIRM
if [[ "$START_CONFIRM" != "yes" ]]; then
  echo "Script terminated."
  exit 0
fi

# ===================== Input Directory and BAM =====================
echo "Enter the FULL path to the directory containing the BAM file without read groups:"
read INPUT_DIR
if [[ ! -d "$INPUT_DIR" ]]; then
  echo "Error: directory not found."
  exit 1
fi

echo "Files found in the directory:"
ls "$INPUT_DIR"
echo "Enter the exact name of the BAM file to be processed:"
read BAM_FILE

# ===================== Read Group Info =====================
echo "Now enter the information required to properly identify the sample in GATK."
echo "These values will be used in the BAM file header."
echo ""
echo "RGID: Read Group ID. This could be the file or sample name."
read -p "Enter RGID (e.g., patient01): " RGID

echo "RGLB: Library identifier. If unsure, use something simple like 'lib1'."
read -p "Enter RGLB (e.g., lib1): " RGLB

echo "RGPL: Sequencing platform used. Choose one of the following options:"
echo "1 = ILLUMINA"
echo "2 = IONTORRENT"
echo "3 = PACBIO"
echo "4 = ONT"
echo "5 = BGI"
read -p "Enter the corresponding number: " RGPL_CHOICE
case $RGPL_CHOICE in
  1) RGPL="ILLUMINA";;
  2) RGPL="IONTORRENT";;
  3) RGPL="PACBIO";;
  4) RGPL="ONT";;
  5) RGPL="BGI";;
  *) echo "Invalid option."; exit 1;;
esac

echo "RGPU: Platform unit. You can make something up, like 'unit1'."
read -p "Enter RGPU: " RGPU

echo "RGSM: Sample name. Use the same name you will use throughout the pipeline."
read -p "Enter RGSM (e.g., patient01): " RGSM

# ===================== Docker User =====================
echo "Do you want to run Docker with a specific user (--user UID)?"
echo "Enter the UID (e.g., 1006), or 'no' to run as root."
read CUSTOM_UID
if [[ "$CUSTOM_UID" =~ ^[0-9]+$ ]]; then
  USE_USER="--user $CUSTOM_UID"
else
  USE_USER=""
fi

# ===================== Resource Allocation =====================
echo "=============================================================="
echo "Note: This tool uses minimal CPU and RAM, but it is best to allocate safely."
echo "Avoid using 100% of system resources to ensure stability."
echo "=============================================================="
TOTAL_THREADS=$(nproc)
TOTAL_RAM=$(free -g | awk '/Mem:/ {print $2}')
echo "Your system has $TOTAL_THREADS threads and ${TOTAL_RAM}GB of RAM available."
echo "How many threads would you like to use? (default = 1)"
read NUM_THREADS
echo "How many GB of RAM would you like to allocate? (default = 2GB)"
echo "Note: Your system has ${TOTAL_RAM}GB. A good choice might be $((TOTAL_RAM - 2))GB."
read MEMORY_RAM

# ===================== Output Directory =====================
TIMESTAMP=$(date +"%d-%m-%Y_%Hh%Mm")
OUTPUT_DIR="$INPUT_DIR/5-output_add_readgroups_$TIMESTAMP"
mkdir -p "$OUTPUT_DIR"
OUTPUT_BAM="$OUTPUT_DIR/${BAM_FILE%.bam}_with_RG.bam"
LOG_FILE="$OUTPUT_DIR/${BAM_FILE%.bam}_addRG.log"

# ===================== Final Summary =====================
echo "================ FINAL SUMMARY ================"
echo "Input BAM file: $BAM_FILE"
echo "Input directory: $INPUT_DIR"
echo "RGID: $RGID"
echo "RGLB: $RGLB"
echo "RGPL: $RGPL"
echo "RGPU: $RGPU"
echo "RGSM: $RGSM"
echo "Docker user: ${CUSTOM_UID:-root}"
echo "Threads: $NUM_THREADS"
echo "RAM: ${MEMORY_RAM}GB"
echo "Output BAM: $OUTPUT_BAM"
echo "Log file: $LOG_FILE"
echo "==============================================="
echo "Ready to start? (yes/no)"
read FINAL_CONFIRM
if [[ "$FINAL_CONFIRM" != "yes" ]]; then
  echo "Execution cancelled."
  exit 0
fi

# ===================== Execution =====================
echo "Starting Read Group insertion..."
START_TIME=$(date +%s)

docker run --rm $USE_USER \
  -v "$INPUT_DIR":/data \
  -v "$OUTPUT_DIR":/output \
  -e _JAVA_OPTIONS="-Xmx${MEMORY_RAM}g" \
  --cpus=$NUM_THREADS \
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
echo "Process completed! BAM with read groups saved at:"
echo "$OUTPUT_BAM"
echo "Total time: ${HOURS}h ${MINUTES}m ${SECONDS}s"
echo "=============================================================="
exit 0
