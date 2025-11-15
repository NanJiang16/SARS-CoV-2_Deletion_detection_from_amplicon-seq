#!/usr/bin/env bash
set -euo pipefail

#######################
# User-configurable variables
#######################

# Sample name (used as base for filenames)
SAMPLE_NAME="sample1"

# Input FASTQ files
R1_FASTQ="${SAMPLE_NAME}_R1.fastq.gz"
R2_FASTQ="${SAMPLE_NAME}_R2.fastq.gz"

# Paths to tools and reference files
TRIMMOMATIC_JAR="/path/to/Trimmomatic-0.39/trimmomatic-0.39.jar"
TRIMMOMATIC_ADAPTERS="/path/to/TruSeq3-PE.fa"

REF_FASTA="SARS-COV2_Reference_MN908947.3_padded.fasta"
PRIMER_BED="SARs-CoV-2_v5.3.2_400.primer.bed"

# ViReMa / pipeline parameters
SEED_LENGTH=20
MICROINDEL_LENGTH=5
MINLEN_TRIMMOMATIC=75
MIN_DELETION_LENGTH=6
MAX_OVERHANG_PRIMER_FRAC=1
MIN_ALIGNED_LENGTH=75

# Output naming
OUTDIR="ViReMa25_SARS2_${SAMPLE_NAME}"
OUTPUT_TAG="ViReMa25_SARS2_${SAMPLE_NAME}"

#######################
# Script begins
#######################

echo ">>> Running pipeline for sample: ${SAMPLE_NAME}"
mkdir -p "${OUTDIR}"

########################################
# 1. Adapter/quality trimming (Trimmomatic)
########################################
echo ">>> Step 1: Trimming adapters and low-quality bases with Trimmomatic"

java -jar "${TRIMMOMATIC_JAR}" PE \
  "${R1_FASTQ}" "${R2_FASTQ}" \
  "output_${SAMPLE_NAME}_paired_R1.fastq" "output_${SAMPLE_NAME}_unpaired_R1.fastq" \
  "output_${SAMPLE_NAME}_paired_R2.fastq" "output_${SAMPLE_NAME}_unpaired_R2.fastq" \
  ILLUMINACLIP:"${TRIMMOMATIC_ADAPTERS}":2:30:10:2:True \
  LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:${MINLEN_TRIMMOMATIC}

########################################
# 2. Interleave paired-end reads
########################################
echo ">>> Step 2: Interleaving paired-end reads"

python3 interleave_paired_end_fastq.py \
  "output_${SAMPLE_NAME}_paired_R1.fastq" \
  "output_${SAMPLE_NAME}_paired_R2.fastq" \
  > "interleaved_${SAMPLE_NAME}.fastq"

########################################
# 3. Run ViReMa
########################################
echo ">>> Step 3: Running ViReMa"

python3 ViReMa.py \
  "${REF_FASTA}" \
  "interleaved_${SAMPLE_NAME}.fastq" \
  "${OUTDIR}/${OUTPUT_TAG}_recombinations.sam" \
  --Output_Dir "${OUTDIR}" \
  --Output_Tag "${OUTPUT_TAG}" \
  --Seed "${SEED_LENGTH}" \
  -BED \
  --MicroInDel_Length "${MICROINDEL_LENGTH}"

########################################
# 4. Standardize alignments
########################################
echo ">>> Step 4: Standardizing alignments"

python3 standardize_alignments.py \
  "${REF_FASTA}" \
  "${OUTDIR}/${OUTPUT_TAG}_recombinations.sam" \
  > "${OUTDIR}/${OUTPUT_TAG}_standardized_recombinations.bam"

########################################
# 5. Filter alignments
########################################
echo ">>> Step 5: Filtering alignments"

python3 filter_aligned_reads.py \
  "${PRIMER_BED}" \
  "${OUTDIR}/${OUTPUT_TAG}_standardized_recombinations.bam" \
  --min-deletion-length "${MIN_DELETION_LENGTH}" \
  --max-overhang-primer-frac "${MAX_OVERHANG_PRIMER_FRAC}" \
  --min-aligned-length "${MIN_ALIGNED_LENGTH}" \
  --virema \
  > "${OUTDIR}/filtered_${OUTPUT_TAG}_standardized_recombinations.bam"

########################################
# 6. Sort filtered BAM (for deletion extraction)
########################################
echo ">>> Step 6: Sorting filtered BAM"

samtools sort \
  "${OUTDIR}/filtered_${OUTPUT_TAG}_standardized_recombinations.bam" \
  -o "${OUTDIR}/filtered_${OUTPUT_TAG}_standardized_recombinations.sorted.bam"

########################################
# 7. Extract deletions
########################################
echo ">>> Step 7: Extracting deletions"

python3 extract_deletions.py \
  "${OUTDIR}/filtered_${OUTPUT_TAG}_standardized_recombinations.sorted.bam" \
  --primer-bed "${PRIMER_BED}" \
  --min-deletion-length "${MIN_DELETION_LENGTH}" \
  --virema \
  > "${OUTDIR}/filtered_${OUTPUT_TAG}_standardized_recombinations.sorted.txt"

########################################
# 8. Coverage calculation from original recombinations
########################################
echo ">>> Step 8: Computing coverage from original recombinations"

samtools view -S -b \
  "${OUTDIR}/${OUTPUT_TAG}_recombinations.sam" \
  > "${OUTDIR}/${OUTPUT_TAG}_recombinations.bam"

samtools sort \
  "${OUTDIR}/${OUTPUT_TAG}_recombinations.bam" \
  -o "${OUTDIR}/${OUTPUT_TAG}_recombinations.sorted.bam"

samtools depth -a -m 0 \
  "${OUTDIR}/${OUTPUT_TAG}_recombinations.sorted.bam" \
  > "${OUTDIR}/${OUTPUT_TAG}_recombinations.coverage"

echo ">>> Pipeline finished for sample: ${SAMPLE_NAME}"
echo "    - Deletions table: ${OUTDIR}/filtered_${OUTPUT_TAG}_standardized_recombinations.sorted.txt"
echo "    - Coverage file:  ${OUTDIR}/${OUTPUT_TAG}_recombinations.coverage"
echo ">>> Use deletion_summarization_and_further_filtration.R for final summarization and plotting."
