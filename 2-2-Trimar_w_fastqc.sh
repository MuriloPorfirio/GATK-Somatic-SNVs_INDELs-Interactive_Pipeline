#!/bin/bash

echo "Script started."
echo "If you pasted this script into the terminal, press Enter to continue (this prevents accidental auto-execution)."
read

echo "How many patients do you want to process? (or type 'exit' to quit)"
read N
if [[ "$N" == "exit" ]]; then exec bash; fi
if ! [[ "$N" =~ ^[0-9]+$ ]]; then
  echo "Invalid value. Please enter an integer."
  exec bash
fi

TOTAL_FILES=$((N * 2))
echo "This corresponds to $TOTAL_FILES files (R1 and R2 per patient). Is that correct? (yes/no or 'exit')"
read CONFIRM_TOTAL
if [[ "$CONFIRM_TOTAL" == "exit" ]]; then exec bash; fi
if [[ "$CONFIRM_TOTAL" != "yes" ]]; then
  echo "Script aborted. Please check your input before proceeding."
  exec bash
fi

echo "Are all files in the same directory? (yes/no or 'exit')"
read SAME_DIR
if [[ "$SAME_DIR" == "exit" ]]; then exec bash; fi
if [[ "$SAME_DIR" != "yes" ]]; then
  echo "Please consolidate all files into a single directory before continuing."
  exec bash
fi

echo "Enter the full path to the directory with the files (no trailing slash, or type 'exit' to quit):"
read INPUT_DIR
if [[ "$INPUT_DIR" == "exit" ]]; then exec bash; fi
if [[ ! -d "$INPUT_DIR" ]]; then
  echo "Error: directory '$INPUT_DIR' not found."
  exec bash
fi

echo "Listing .fastq or .fastq.gz files in '$INPUT_DIR'..."
find "$INPUT_DIR" -maxdepth 1 -type f \( -name "*.fastq" -o -name "*.fastq.gz" \)

echo "Do you want to process all the files listed above? (yes/no or 'exit')"
read PROCESS_ALL
if [[ "$PROCESS_ALL" == "exit" ]]; then exec bash; fi

if [[ "$PROCESS_ALL" != "yes" ]]; then
  echo "Enter the EXACT names of the files you want to trim, separated by space (or 'exit' to quit):"
  read -a FILES
  if [[ "${FILES[0]}" == "exit" ]]; then exec bash; fi
else
  FILES=($(find "$INPUT_DIR" -maxdepth 1 -type f \( -name "*.fastq" -o -name "*.fastq.gz" \) -printf "%f\n"))
fi

# Create output folder
TIMESTAMP=$(date +"%d-%m-%Y_%Hh%Mm")
OUTPUT_DIR="$INPUT_DIR/2-trimmed_files_$TIMESTAMP"
mkdir -p "$OUTPUT_DIR"

echo "Trimmed files will be saved in: $OUTPUT_DIR"
echo "Do you want to continue? (yes/no or 'exit')"
read CONTINUE
if [[ "$CONTINUE" == "exit" ]]; then exec bash; fi
if [[ "$CONTINUE" != "yes" ]]; then
  echo "Process cancelled."
  exec bash
fi

# Detect system resources
TOTAL_CPUS=$(nproc)
SAFE_CPUS=$((TOTAL_CPUS / 2))
TOTAL_RAM_GB=$(free -g | awk '/^Mem:/ {print $2}')
DEFAULT_RAM_GB=4
SUGGESTED_RAM_GB=$((TOTAL_RAM_GB / 2))

echo ""
echo "This server has $TOTAL_CPUS threads available."
echo "It is recommended to use up to $SAFE_CPUS threads to avoid impacting other users."
echo "How many threads do you want to use? (or 'exit')"
read THREADS
if [[ "$THREADS" == "exit" ]]; then exec bash; fi
if ! [[ "$THREADS" =~ ^[0-9]+$ ]]; then
  echo "Invalid number of threads."
  exec bash
fi

echo ""
echo "This script typically uses around ${DEFAULT_RAM_GB}GB of RAM."
echo "The server has $TOTAL_RAM_GB GB. We suggest using up to ${SUGGESTED_RAM_GB}GB."
echo "How much RAM would you like to allocate (in GB)? (or 'exit')"
read RAM
if [[ "$RAM" == "exit" ]]; then exec bash; fi
if ! [[ "$RAM" =~ ^[0-9]+$ ]]; then
  echo "Invalid value for RAM."
  exec bash
fi

# Docker --user
echo ""
echo "Do you want to run Docker with a specific user (--user UID)?"
echo "Enter the desired UID (e.g., 1006), or 'no' to continue as default:"
read CUSTOM_UID
if [[ "$CUSTOM_UID" == "exit" ]]; then exec bash; fi
if [[ "$CUSTOM_UID" =~ ^[0-9]+$ ]]; then
  USE_USER="--user $CUSTOM_UID"
  echo "Docker will run with: $USE_USER"
elif [[ "$CUSTOM_UID" == "no" ]]; then
  USE_USER=""
  echo "Docker will run as default user (root inside the container)."
else
  echo "Invalid input. Please enter a UID number or 'no'."
  exec bash
fi

# Final summary
echo ""
echo "========= FINAL SUMMARY ========="
echo "Input directory: $INPUT_DIR"
echo "Output directory: $OUTPUT_DIR"
echo "Files to be processed:"
for file in "${FILES[@]}"; do
  echo "- $file"
done
echo ""
echo "Threads: $THREADS"
echo "RAM (informational): $RAM GB"
echo "Docker UID: ${CUSTOM_UID}"
echo "================================="
echo "Shall we begin? (yes/no or 'exit')"
read FINAL_CONFIRM
if [[ "$FINAL_CONFIRM" == "exit" ]]; then exec bash; fi
if [[ "$FINAL_CONFIRM" != "yes" ]]; then
  echo "Execution cancelled."
  exec bash
fi

# Start trimming
echo "Starting trimming with Docker..."

i=0
while [[ $i -lt ${#FILES[@]} ]]; do
  R1="${FILES[$i]}"
  R2="${FILES[$((i+1))]}"
  echo "Processing pair:"
  echo "R1: $R1"
  echo "R2: $R2"

  docker run --rm \
    -v "$INPUT_DIR":/data \
    $USE_USER \
    biowardrobe2/trimgalore:v0.4.4 \
    trim_galore --paired "/data/$R1" "/data/$R2" -o "/data/2-trimmed_files_$TIMESTAMP"

  if [[ $? -ne 0 ]]; then
    echo "Error processing: $R1 and $R2"
  else
    echo "Completed: $R1 and $R2"
  fi

  i=$((i + 2))
done

echo "------------------------------------------------------"
echo "Trimming process finished."
echo "Results are in: $OUTPUT_DIR"
echo "------------------------------------------------------"

# Optionally run FastQC
echo ""
echo "Would you like to run FastQC on the trimmed files now? (yes/no or 'exit')"
read QC_CONFIRM
if [[ "$QC_CONFIRM" == "exit" ]]; then exec bash; fi

if [[ "$QC_CONFIRM" == "yes" ]]; then
  echo "Running FastQC on trimmed files..."
  docker run --rm \
    -v "$OUTPUT_DIR":/data \
    biocontainers/fastqc:v0.11.9_cv8 \
    bash -c "fastqc -o /data /data/*_val_1.fq.gz /data/*_val_2.fq.gz"

  if [[ $? -ne 0 ]]; then
    echo "There was an error running FastQC."
  else
    echo "FastQC completed successfully. Reports are in the output directory."
  fi
else
  echo "FastQC skipped as requested."
fi

echo "Script finished."
exec bash
