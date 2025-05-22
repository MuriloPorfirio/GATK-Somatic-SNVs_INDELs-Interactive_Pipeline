#!/bin/bash

# Start
echo "========================================="
echo "         Script to Run MultiQC           "
echo "========================================="
echo ""
echo "If you pasted this script into the terminal, press Enter to continue."
read

# Ask if Docker should use a specific UID
echo "Do you want to run Docker with a specific user (--user UID)?"
echo "Enter the UID number (e.g., 1006), or 'no' to run as root, or 'exit' to quit."
read CUSTOM_UID
if [[ "$CUSTOM_UID" == "exit" ]]; then exit 0; fi

if [[ "$CUSTOM_UID" =~ ^[0-9]+$ ]]; then
  USE_USER="--user $CUSTOM_UID"
  echo "Docker will run with: $USE_USER"
elif [[ "$CUSTOM_UID" == "no" ]]; then
  USE_USER=""
  echo "Docker will run as root inside the container."
else
  echo "Invalid input. Exiting."
  exit 1
fi

# Input directory
echo "Now enter the FULL path to the directory where the files are located (all files must be in the same directory)."
echo "To be sure, navigate to the folder and type 'pwd', then copy and paste the path here."
read INPUT_DIR
if [[ "$INPUT_DIR" == "exit" ]]; then exit 0; fi

if [[ ! -d "$INPUT_DIR" ]]; then
  echo "Error: Directory '$INPUT_DIR' not found."
  exit 1
fi

# Show available files
echo "Listing all files in '$INPUT_DIR'..."
FILES_AVAILABLE=($(ls "$INPUT_DIR"))
for file in "${FILES_AVAILABLE[@]}"; do
  echo "- $file"
done

# Get the files to be used by MultiQC
echo "Now copy and paste the names of the files you want to include in MultiQC, separated by space."
echo "NOTE: Only .zip files are allowed (example below):"
echo "patient1_R1_fastqc.zip patient1_R2_fastqc.zip patient2_R1_fastqc.zip patient2_R2_fastqc.zip"
echo "(Use exactly the file names as shown above, including case and extensions!)"
read -a FILES_SELECTED

# Basic validation
if [[ -z "${FILES_SELECTED[*]}" ]]; then
  echo "No files provided. Exiting."
  exit 1
fi

# Check file validity
for file in "${FILES_SELECTED[@]}"; do
  if [[ "$file" == *.html ]]; then
    echo "Error: '$file' is an .html file. MultiQC requires .zip files, not .html."
    echo "Exiting."
    exit 1
  fi
  if [[ "$file" != *.zip ]]; then
    echo "Error: '$file' is not a .zip file. Only .zip files are accepted."
    echo "Exiting."
    exit 1
  fi
  if [[ ! -f "$INPUT_DIR/$file" ]]; then
    echo "Error: File '$file' not found in the specified directory."
    echo "Exiting."
    exit 1
  fi

  # Add full path for container
  FILES_TO_ANALYZE+=("/data/$file")
done

# Create output directory
TIMESTAMP=$(date +"%d-%m-%Y_%Hh%Mm")
OUTPUT_DIR="$INPUT_DIR/3-multiqc_output_$TIMESTAMP"
mkdir -p "$OUTPUT_DIR"

# Final summary
echo "========= FINAL SUMMARY ========="
echo "Input directory: $INPUT_DIR"
echo "Files selected for MultiQC:"
for file in "${FILES_SELECTED[@]}"; do
  echo "- $file"
done
echo ""
echo "Output directory: $OUTPUT_DIR"
echo "Docker UID: ${CUSTOM_UID}"
echo "================================="
echo ""
echo "Is everything correct? (yes/no)"
read FINAL_CONFIRM
if [[ "$FINAL_CONFIRM" == "exit" ]]; then exit 0; fi
if [[ "$FINAL_CONFIRM" != "yes" ]]; then
  echo "Execution cancelled."
  exit 0
fi

# Run MultiQC
echo "Starting MultiQC execution..."
docker run --rm \
  -v "$INPUT_DIR":/data \
  -v "$OUTPUT_DIR":/output \
  $USE_USER \
  ewels/multiqc \
  multiqc /data --outdir /output --force

# Done
echo ""
echo "==================================================="
echo "MultiQC completed."
echo "Report available in: $OUTPUT_DIR"
echo "==================================================="
exit 0
