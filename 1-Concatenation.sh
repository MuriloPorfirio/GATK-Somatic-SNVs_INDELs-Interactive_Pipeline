#!/bin/bash

echo "Starting FASTQ file concatenation script."
echo "Press Enter to confirm you intentionally pasted this script into the terminal."
read

# Initial question
echo "How many patients will you be working with? (or type 'exit' to quit)"
read PATIENTS
if [[ "$PATIENTS" == "exit" ]]; then exec bash; fi
if ! [[ "$PATIENTS" =~ ^[0-9]+$ ]]; then
  echo "Invalid input. Please enter an integer number."
  exec bash
fi

# Paired files?
echo "Are the files for these patients paired-end? (yes/no or 'exit')"
read PAIRED
if [[ "$PAIRED" == "exit" ]]; then exec bash; fi

if [[ "$PAIRED" != "yes" ]]; then
  echo "This script only works with paired-end files. Exiting."
  exec bash
fi

# Automatic calculation
TOTAL_FILES=$((PATIENTS * 4))
R1_FILES_COUNT=$((PATIENTS * 2))
R2_FILES_COUNT=$((PATIENTS * 2))

echo ""
echo "There will be $TOTAL_FILES total files:"
echo "- $R1_FILES_COUNT files for R1"
echo "- $R2_FILES_COUNT files for R2"
echo ""

# File path
echo "Enter the FULL path to the directory containing the files (no trailing slash, or 'exit' to quit):"
read FILE_PATH
if [[ "$FILE_PATH" == "exit" ]]; then exec bash; fi
if [[ ! -d "$FILE_PATH" ]]; then
  echo "Directory '$FILE_PATH' not found."
  exec bash
fi

echo "Files found in $FILE_PATH:"
ls "$FILE_PATH"
echo ""
echo "All files must be in this directory."
echo "If not, type 'exit' and move the files to a single location."
echo "Do you want to continue? (yes/no or 'exit')"
read CONTINUE
if [[ "$CONTINUE" == "exit" ]]; then exec bash; fi
if [[ "$CONTINUE" != "yes" ]]; then echo "Cancelled."; exec bash; fi

# R1 files
echo "Enter the EXACT names of the $R1_FILES_COUNT R1 files, separated by space:"
read -a R1_FILES
if [[ "${#R1_FILES[@]}" -ne "$R1_FILES_COUNT" ]]; then
  echo "Incorrect number of R1 files provided."
  exec bash
fi

# R2 files
echo "Enter the EXACT names of the $R2_FILES_COUNT R2 files, separated by space:"
read -a R2_FILES
if [[ "${#R2_FILES[@]}" -ne "$R2_FILES_COUNT" ]]; then
  echo "Incorrect number of R2 files provided."
  exec bash
fi

# Output folder
TIMESTAMP=$(date +"%d-%m-%Y_%Hh%Mm")
OUTPUT_DIR="$FILE_PATH/1-concatenated_fastq_$TIMESTAMP"
mkdir -p "$OUTPUT_DIR"

echo ""
echo "The concatenated files will be saved in: $OUTPUT_DIR"
echo "Do you want to proceed with concatenation? (yes/no or 'exit')"
read CONFIRM_FINAL
if [[ "$CONFIRM_FINAL" == "exit" ]]; then exec bash; fi
if [[ "$CONFIRM_FINAL" != "yes" ]]; then echo "Cancelled."; exec bash; fi

# Concatenation
echo ""
echo "Starting concatenation..."
for ((i=0; i<PATIENTS; i++)); do
  R1_1=${R1_FILES[$((i*2))]}
  R1_2=${R1_FILES[$((i*2+1))]}
  R2_1=${R2_FILES[$((i*2))]}
  R2_2=${R2_FILES[$((i*2+1))]}

  PREFIX=$(echo "${R1_1}" | cut -d'_' -f1-2)

  echo "Concatenating files for sample: $PREFIX"
  cat "$FILE_PATH/$R1_1" "$FILE_PATH/$R1_2" > "$OUTPUT_DIR/${PREFIX}_R1_combined.fastq.gz"
  cat "$FILE_PATH/$R2_1" "$FILE_PATH/$R2_2" > "$OUTPUT_DIR/${PREFIX}_R2_combined.fastq.gz"

  echo "Concatenation completed: $PREFIX"
done

echo ""
echo "All files were successfully concatenated."
echo "Final files are located in: $OUTPUT_DIR"
exec bash
