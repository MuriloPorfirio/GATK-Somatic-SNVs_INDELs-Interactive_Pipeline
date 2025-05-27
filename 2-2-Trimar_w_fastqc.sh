# Configure advanced Trim Galore parameters
echo ""
echo "Set the base quality cutoff for trimming (default is 20)."
echo "This trims low-quality ends from reads before adapter removal."
echo "Typical values range from 15 to 30."
echo "Recommended: 20 for WES; 25 for WGS with somatic calls."
echo "Enter your desired quality value or press Enter to use default:"
read QUALITY
if [[ -z "$QUALITY" ]]; then
  QUALITY=20
fi

echo ""
echo "Set the minimum read length to keep after trimming (default is 20)."
echo "Shorter reads are discarded. Acceptable values: 20–75."
echo "Recommended: 30 for WES; 50 for WGS with somatic calls."
echo "Enter your desired minimum length or press Enter to use default:"
read MIN_LENGTH
if [[ -z "$MIN_LENGTH" ]]; then
  MIN_LENGTH=20
fi

echo ""
echo "Would you like to manually specify the adapter type? (illumina/nextera/auto)"
echo "Trim Galore detects adapters automatically, so you can safely choose 'auto'."
echo "If you choose 'illumina' or 'nextera', this will force that adapter trimming mode."
read ADAPTER_TYPE
if [[ "$ADAPTER_TYPE" == "illumina" ]]; then
  ADAPTER_OPTION="--illumina"
elif [[ "$ADAPTER_TYPE" == "nextera" ]]; then
  ADAPTER_OPTION="--nextera"
else
  ADAPTER_OPTION=""
fi

echo ""
echo "Do you want to run FastQC after trimming? (yes/no)"
read RUN_FASTQC
if [[ "$RUN_FASTQC" == "yes" ]]; then
  FASTQC_OPTION="--fastqc"
else
  FASTQC_OPTION="--no_fastqc"
fi

# Recap of trimming parameters
echo ""
echo "======= TRIM PARAMETERS SUMMARY ======="
echo "Quality cutoff: $QUALITY"
echo "Minimum read length: $MIN_LENGTH"
echo "Adapter type: ${ADAPTER_TYPE:-auto-detect}"
echo "Run FastQC: $RUN_FASTQC"
echo "======================================="

# Start trimming with time tracking
START_TIME=$(date +%s)
LOG_FILE="$OUTPUT_DIR/trimming_log_$(date +"%Y-%m-%d_%H%M").log"
echo "Log file: $LOG_FILE"
echo "Starting trimming..." | tee -a "$LOG_FILE"

i=0
while [[ $i -lt ${#FILES[@]} ]]; do
  R1="${FILES[$i]}"
  R2="${FILES[$((i+1))]}"
  echo "Processing: $R1 and $R2" | tee -a "$LOG_FILE"

  docker run --rm \
    -v "$INPUT_DIR":/data \
    $USE_USER \
    biowardrobe2/trimgalore:v0.4.4 \
    trim_galore --paired \
      --quality "$QUALITY" \
      --length "$MIN_LENGTH" \
      $ADAPTER_OPTION \
      $FASTQC_OPTION \
      "/data/$R1" "/data/$R2" \
      -o "/data/2-trimmed_files_$TIMESTAMP" 2>&1 | tee -a "$LOG_FILE"

  if [[ ${PIPESTATUS[0]} -ne 0 ]]; then
    echo "Error processing $R1 and $R2" | tee -a "$LOG_FILE"
  else
    echo "Completed $R1 and $R2" | tee -a "$LOG_FILE"
  fi

  i=$((i + 2))
done

END_TIME=$(date +%s)
TOTAL_TIME=$((END_TIME - START_TIME))
echo "Trimming completed." | tee -a "$LOG_FILE"
echo "Total runtime: $TOTAL_TIME seconds" | tee -a "$LOG_FILE"

echo "------------------------------------------------------"
echo "All trimming steps are complete."
echo "Results saved in: $OUTPUT_DIR"
echo "Full log: $LOG_FILE"
echo "Total processing time: $TOTAL_TIME seconds"
echo "------------------------------------------------------"
