#!/bin/bash

# ===================== Introduction =====================
echo "=============================================================="
echo "        Interactive Script for GATK MarkDuplicates"
echo "=============================================================="
echo ""
echo "GATK MarkDuplicates identifies and marks PCR duplicate reads"
echo "to avoid biases in variant analysis."
echo ""
echo "IMPORTANCE OF THIS TOOL:"
echo " - Removes biases caused by PCR duplicates."
echo " - Improves the accuracy of variant calling."
echo " - Maintains real (not artificially inflated) read depth."
echo ""
echo "REQUIREMENTS BEFORE RUNNING THIS SCRIPT:"
echo " - BAM files must be aligned and sorted (produced after BWA alignment and Samtools sorting)."
echo ""
echo "Shall we continue? (yes/no)"
read START_CONFIRM
if [[ "$START_CONFIRM" != "yes" ]]; then
  echo "Script terminated."
  exit 0
fi

# ===================== Input BAM Files =====================
echo "Enter the FULL path to the directory containing the aligned BAM files:"
echo "(Use 'pwd' to copy the correct path, right-click to paste.)"
read INPUT_DIR
if [[ "$INPUT_DIR" == "exit" ]]; then exit 0; fi
if [[ ! -d "$INPUT_DIR" ]]; then
  echo "Error: Directory '$INPUT_DIR' not found."
  exit 1
fi

echo "Files found in the selected directory:"
ls "$INPUT_DIR"
echo ""
echo "Enter the exact names of the BAM files you want to process, separated by space."
echo "Example: patient1.bam patient2.bam patient3.bam"
read -a BAM_FILES
if [[ "$BAM_FILES" == "exit" ]]; then exit 0; fi

# ===================== Docker User =====================
echo "Do you want to run Docker with a specific user (--user UID)?"
echo "Enter the UID (e.g.: 1006), or 'no' to run as root, or 'exit' to quit:"
read CUSTOM_UID
if [[ "$CUSTOM_UID" == "exit" ]]; then exit 0; fi
if [[ "$CUSTOM_UID" =~ ^[0-9]+$ ]]; then
  USE_USER="--user $CUSTOM_UID"
else
  USE_USER=""
fi

# ===================== System Resources =====================
echo "GATK MarkDuplicates uses CPU. The default is to use 1 thread."
TOTAL_THREADS=$(nproc)
echo "Your system has $TOTAL_THREADS available threads."
echo "How many threads would you like to use? (recommended: leave some free, e.g.: $((TOTAL_THREADS - 2)))"
read NUM_THREADS

DEFAULT_RAM=8
echo "GATK MarkDuplicates typically uses around ${DEFAULT_RAM}GB of RAM."
TOTAL_RAM=$(free -g | awk '/Mem:/ {print $2}')
echo "Your system has ${TOTAL_RAM}GB of RAM."
echo "How many GB of RAM would you like to allocate? (recommended: $((TOTAL_RAM - 2)))"
read MEMORY_RAM

# ===================== Create Output Directory =====================
TIMESTAMP=$(date +"%d-%m-%Y_%Hh%Mm")
OUTPUT_DIR="$INPUT_DIR/6-output_marked_duplicates_$TIMESTAMP"
mkdir -p "$OUTPUT_DIR"

# ===================== Summary =====================
echo "================ SUMMARY ================"
echo "Selected files:"
for file in "${BAM_FILES[@]}"; do echo "- $file"; done
echo ""
echo "BAM directory: $INPUT_DIR"
echo "Docker user: ${CUSTOM_UID:-root}"
echo "Threads to use: $NUM_THREADS"
echo "RAM to use: ${MEMORY_RAM}GB"
echo "Output will be saved in: $OUTPUT_DIR"
echo "========================================="
echo "Shall we begin? (yes/no)"
read FINAL_CONFIRM
if [[ "$FINAL_CONFIRM" != "yes" ]]; then
  echo "Execution cancelled."
  exit 0
fi

# ===================== Execution =====================
echo "Starting GATK MarkDuplicates execution..."
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
  echo "Time to process ${SAMPLE_NAME}: ${SAMPLE_DURATION} seconds." | tee -a "$LOG_FILE"
done

SCRIPT_END_TIME=$(date +%s)
TOTAL_DURATION=$((SCRIPT_END_TIME - SCRIPT_START_TIME))
TOTAL_HOURS=$((TOTAL_DURATION / 3600))
TOTAL_MINUTES=$(((TOTAL_DURATION % 3600) / 60))
TOTAL_SECONDS=$((TOTAL_DURATION % 60))

echo ""
echo "=============================================================="
echo "GATK MarkDuplicates completed! BAM files and logs are located at:"
echo "$OUTPUT_DIR"
echo "Total execution time: ${TOTAL_HOURS}h ${TOTAL_MINUTES}m ${TOTAL_SECONDS}s."
echo "=============================================================="
exit 0
