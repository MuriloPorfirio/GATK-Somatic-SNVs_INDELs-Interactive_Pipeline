#  Interactive Pipeline for SNVs and INDELs – Processing Illumina NGS FASTQs

This repository contains Bash scripts for preprocessing genomic sequencing data (FASTQ files).

---

## ⚠️ Important Guidelines Before Use

## Docker-Based Workflow

Most steps in this pipeline are **containerized using Docker**. This ensures reproducibility and avoids software conflicts.  
Before running the scripts, make sure the following Docker images are installed locally (or available for pull):

biowardrobe2/trimgalore:v0.4.4
biocontainers/fastqc:v0.11.9_cv8
ewels/multiqc
biocontainers/bwa:v0.7.17_cv1
staphb/samtools:latest
broadinstitute/gatk:latest

### Interactive and long scripts

All scripts in this repository are **interactive**, prompting the user with several questions. Because of that, they are visually long and **should not be pasted directly into a terminal within a screen session**, because:

- The content might be **executed instantly**, without time to respond
- This can freeze or crash the session without warning

> **Recommended usage:**
> 1. Create a `.sh` file inside your `screen` session using `nano`:
>    ```bash
>    nano my_script.sh
>    ```
> 2. Paste the script content into the nano editor
> 3. Save it (`Ctrl+O`, `Enter`, `Ctrl+X`)
> 4. Make the script executable:
>    ```bash
>    chmod +x my_script.sh
>    ```
> 5. Run it inside the screen with:
>    ```bash
>    ./my_script.sh
>    ```

---

### Correct order: R1 must always come before R2

**Throughout all steps of the pipeline up to the final alignment**, it is crucial to ensure that **R1 files are always listed before R2**, especially in commands such as:

- Trimming (`trim_galore`)
- FastQC/MultiQC
- Alignment (Bowtie2, Parabricks, etc.)

> An incorrect order may affect analysis quality or produce silent errors.

## Best Practices

- Name your files clearly (e.g., `_R1_`, `_R2_`, `_val_1`, `_val_2`)
- Run the scripts **in the proper order**
- Use `screen` or `tmux` with caution
- Always review the generated `.html` reports
- Always maintain the correct read order: **R1 first, R2 second**

---

## Author
Murilo Porfírio de Aguiar and GPT  
murilo.porfirio@yahoo.com
