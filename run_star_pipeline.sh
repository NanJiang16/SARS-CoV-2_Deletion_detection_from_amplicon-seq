#!/usr/bin/env bash
set -euo pipefail

#######################
# User-configurable variables
#######################

# Sample identifiers
SRA="SRRXXXXXXX"          # Folder name and FASTQ prefix
NAME="sample_star_"       # STAR output prefix (will be used by STAR and downstream)

# Input FASTQ files (untrimmed, paired-end)
R1_FASTQ="./${SRA}/${SRA}_1.fastq"
R2_FASTQ="./${SRA}/${SRA}_2.fastq"

# Trim Galore
TRIMGALORE_BIN="$HOME/TrimGalore-0.4.3/trim_galore"
TRIMGALORE_STRINGENCY=3
TRIMGALORE_QSCORE=30
TRIMGALORE_ERROR_RATE=0.10
TRIMGALORE_MINLEN=15

# STAR reference and parameters
GENOME_DIR="./Genome_Dir"                 # STAR genome directory
STAR_REF_FASTA="NC_045512.2.fasta"       # Reference fasta used in standardize_alignments.py

# Primer BED file
PRIMER_BED="ARTIC_primers_v3.bed"

# Filtering parameters
MIN_DELETION_LENGTH=20
MAX_OVERHANG_PRIMER_FRAC=1
MIN_ALIGNED_LENGTH=75
MAX_PRIMER_DIST=1

# R summarization script 
RUN_R_SUMMARY=true    # set to false to skip R step
R_SUMMARY_SCRIPT="./deletion_summarization_and_further_filtration.R"

#######################
# Script begins
#######################

echo ">>> Running Trim Galore + STAR pipeline for sample: ${SRA}"

########################################
# 1. Trim adapters and low-quality bases (Trim Galore)
########################################
echo ">>> Step 1: Trimming adapters and low-quality bases with Trim Galore"

"${TRIMGALORE_BIN}" \
  --stringency "${TRIMGALORE_STRINGENCY}" \
  -q "${TRIMGALORE_QSCORE}" \
  -e "${TRIMGALORE_ERROR_RATE}" \
  --length "${TRIMGALORE_MINLEN}" \
  --paired \
  "${R1_FASTQ}" \
  "${R2_FASTQ}"

# Trim Galore output names (default for paired-end):
#   ${SRA}_1_val_1.fq
#   ${SRA}_2_val_2.fq

TRIMMED_R1="./${SRA}_1_val_1.fq"
TRIMMED_R2="./${SRA}_2_val_2.fq"

########################################
# 2. Align trimmed reads with STAR
########################################
echo ">>> Step 2: Aligning trimmed reads with STAR"

STAR \
  --readFilesIn "${TRIMMED_R1}" "${TRIMMED_R2}" \
  --outFileNamePrefix "${NAME}" \
  --genomeDir "${GENOME_DIR}" \
  --outFilterType BySJout \
  --outFilterMultimapNmax 20 \
  --alignSJoverhangMin 8 \
  --alignSJDBoverhangMin 1 \
  --outSJfilterOverhangMin 12 12 12 12 \
  --outSJfilterCountUniqueMin 1 1 1 1 \
  --outSJfilterCountTotalMin 1 1 1 1 \
  --outSJfilterDistToOtherSJmin 0 0 0 0 \
  --outFilterMismatchNmax 999 \
  --outFilterMismatchNoverReadLmax 0.04 \
  --scoreGapNoncan -4 \
  --scoreGapATAC -4 \
  --chimScoreJunctionNonGTAG 0 \
  --chimOutType Junctions WithinBAM HardClip \
  --alignSJstitchMismatchNmax -1 -1 -1 -1 \
  --alignIntronMin 20 \
  --alignIntronMax 1000000 \
  --alignMatesGapMax 1000000

# STAR main SAM output: ${NAME}Aligned.out.sam

########################################
# 3. Standardize alignments (shift deletions downstream)
########################################
echo ">>> Step 3: Standardizing alignments"

python3 standardize_alignments.py \
  "${STAR_REF_FASTA}" \
  "${NAME}Aligned.out.sam" \
  > "${NAME}_standardizationAligned.out.bam"

########################################
# 4. Annotate alignments with primer distances
########################################
echo ">>> Step 4: Annotating alignments with primer information"

python3 annotate_alignment_with_primers.py \
  "${PRIMER_BED}" \
  "${NAME}_standardizationAligned.out.bam" \
  > "${NAME}_annotated_standardization_Aligned.out.bam"

########################################
# 5. Filter alignments (min length, primer pool matching, etc.)
########################################
echo ">>> Step 5: Filtering alignments"

python3 filter_aligned_reads.py \
  "${PRIMER_BED}" \
  "${NAME}_annotated_standardization_Aligned.out.bam" \
  --min-deletion-length "${MIN_DELETION_LENGTH}" \
  --max-overhang-primer-frac "${MAX_OVERHANG_PRIMER_FRAC}" \
  --min-aligned-length "${MIN_ALIGNED_LENGTH}" \
  --primer-pool-matching \
  --max-primer-dist "${MAX_PRIMER_DIST}" \
  > "filtered_${NAME}_annotated_standardizationAligned.out.bam"

########################################
# 6. Sort filtered BAM
########################################
echo ">>> Step 6: Sorting filtered BAM"

samtools sort \
  "filtered_${NAME}_annotated_standardizationAligned.out.bam" \
  -o "filtered_${NAME}_annotated_standardizationAligned.out.sorted.bam"

########################################
# 7. Extract deletions
########################################
echo ">>> Step 7: Extracting deletions"

python3 extract_deletions.py \
  "filtered_${NAME}_annotated_standardizationAligned.out.sorted.bam" \
  --primer-bed "${PRIMER_BED}" \
  --min-deletion-length "${MIN_DELETION_LENGTH}" \
  --ignore-secondary \
  > "filtered_${NAME}_annotated_standardizationAligned.deletion.sorted.txt"

########################################
# 8. Compute coverage
########################################
echo ">>> Step 8: Computing coverage from filtered alignment"

samtools depth -a -m 0 \
  "filtered_${NAME}_annotated_standardizationAligned.out.sorted.bam" \
  > "filtered_${NAME}_annotated_standardizationAligned.out.sorted.coverage"

echo ">>> STAR pipeline finished for sample: ${SRA}"
echo "    - Deletion table: filtered_${NAME}_annotated_standardizationAligned.deletion.sorted.txt"
echo "    - Coverage file:  filtered_${NAME}_annotated_standardizationAligned.out.sorted.coverage"

########################################
# 9. Run R summarization / further filtration 
########################################
if [ "${RUN_R_SUMMARY}" = true ]; then
  echo ">>> Step 9: Running R summarization script: ${R_SUMMARY_SCRIPT}"
  if [ ! -f "${R_SUMMARY_SCRIPT}" ]; then
    echo "ERROR: R summarization script not found: ${R_SUMMARY_SCRIPT}" >&2
    exit 1
  fi
  Rscript "${R_SUMMARY_SCRIPT}"
  echo ">>> R summarization completed for sample(s) in current directory."
else
  echo ">>> Skipping R summarization step (RUN_R_SUMMARY=false)."
fi
