# Deletion detection in SARS-CoV-2 genomes using multiplex-PCR sequencing from COVID-19 patients: elimination of false positives

## Scripts

* `interleave_paired_end_fastq.py` - Interleaves paired-end FASTQ files, appending _1 and _2 to the read names.  To be used for creating input files for Virema.
* `filter_aligned_reads.py` - Filters aligned reads that are likely to have resulted from primer dimers or similar artifacts.
* `extract_deletions.py` - Creates a table of deletions from an alignment file, with information about the primers that flank the deletion.
* `standardize_alignments.py` - Standardizes the deletion positions within alignments.
* `annotate_alignment_with_primers.py` - Annotates an alignment file with information about the closest upstream/downstream primer to each read as well as the fraction of the read's alignment that is to primer regions.
* `extract_primer_info.py` - Extracts primer information from a primer-annotated alignment file.
* `summmarize_deletions.R` - Summarizes read deletions (output from `extract_deletions.py`) and calculates their frequencies using coverage information.

## Python modules

The above scripts rely on the following Python module files:

* `fastq.py` - Classes for reading and writing FASTQ files.
* `dvg.py` - Utility functions for extracting deletions from alignments and comparing with primer information.

## Usage examples 
Filtration and standardization in combination with ViReMa (v0.25) including preprocess with Trimmomatic (0.39)
```bash
java -jar Trimmomatic-0.39/trimmomatic-0.39.jar PE ${name}_R1.fastq.gz ${name}_R2.fastq.gz output_${name}_paired_R1.fastq output_${name}_unpaired_R1.fastq output_${name}_paired_R2.fastq output_${name}_unpaired_R2.fastq ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:2:True LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:75 
python3 interleave_paired_end_fastq.py output_${name}_paired_R1.fastq output_${name}_paired_R2.fastq > interleaved_${name}.fastq
python3 ViReMa.py SARS-COV2_Reference_MN908947.3_padded.fasta interleaved_${name}.fastq ViReMa25_SARS2_${name}_recombinations.sam --Output_Dir ViReMa25_SARS2_${name} --Output_Tag ViReMa25_SARS2_${name} --Seed 20 -BED --MicroInDel_Length 5
python3 standardize_alignments.py SARS-COV2_Reference_MN908947.3_padded.fasta ./ViReMa25_SARS2_${name}/ViReMa25_SARS2_${name}_recombinations.sam > ./ViReMa25_SARS2_${name}/ViReMa25_SARS2_standardized_${name}_recombinations.bam
python3 filter_aligned_reads.py SARs-CoV-2_v5.3.2_400.primer.bed ./ViReMa25_SARS2_${name}/ViReMa25_SARS2_standardized_${name}_recombinations.bam --min-deletion-length 6 --max-overhang-primer-frac 1 --min-aligned-length 75 --virema > ./ViReMa25_SARS2_${name}/filtered_ViReMa25_SARS2_standardized_${name}_recombinations.bam
python3 extract_deletions.py ./ViReMa25_SARS2_${name}/filtered_ViReMa25_SARS2_standardized_${name}_recombinations.sorted.bam --primer-bed SARs-CoV-2_v5.3.2_400.primer.bed --min-deletion-length 6 --virema > ./ViReMa25_SARS2_${name}/filtered_ViReMa25_SARS2_standardized_${name}_recombinations.sorted.txt
samtools view -S -b ViReMa25_SARS2_${name}_recombinations.sam > ViReMa25_SARS2_${name}_recombinations.bam 
samtools sort ViReMa25_SARS2_${name}_recombinations.bam -o ViReMa25_SARS2_${name}_recombinations.sorted.bam
samtools depth -a -m 0 ViReMa25_SARS2_${name}_recombinations.sorted.bam > ViReMa25_SARS2_${name}_recombinations.coverage
```
Filtration and standardization in combination with STAR (v2.7.3a) including preprocess with TrimGalore (0.4.3)
```bash
~/TrimGalore-0.4.3/trim_galore --stringency 3 -q 30 -e .10 --length 15 --paired ./${SRA}/${SRA}_1.fastq ./${SRA}/${SRA}_2.fastq
STAR --readFilesIn ./${SRA}_1_val_1.fq ./${SRA}_2_val_2.fq --outFileNamePrefix ${name} --genomeDir ./Genome_Dir --outFilterType BySJout --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outSJfilterOverhangMin 12 12 12 12 --outSJfilterCountUniqueMin 1 1 1 1 --outSJfilterCountTotalMin 1 1 1 1 --outSJfilterDistToOtherSJmin 0 0 0 0 --outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.04 --scoreGapNoncan -4 --scoreGapATAC -4 --chimScoreJunctionNonGTAG 0 --chimOutType Junctions WithinBAM HardClip --alignSJstitchMismatchNmax -1 -1 -1 -1 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000
python3 standardize_alignments.py NC_045512.2.fasta ${name}Aligned.out.sam > ${name}_standardizationAligned.out.bam
python3 annotate_alignment_with_primers.py ARTIC_primers_v3.bed ${name}_standardizationAligned.out.bam > ${name}_annotated_standardization_Aligned.out.bam
python3 filter_aligned_reads.py ARTIC_primers_v3.bed ${name}_annotated_standardization_Aligned.out.bam --min-deletion-length 20 --max-overhang-primer-frac 1 --min-aligned-length 75 --primer-pool-matching --max-primer-dist 1 > filtered_${name}_annotated_standardizationAligned.out.bam
samtools sort filtered_${name}_annotated_standardizationAligned.out.bam -o filtered_${name}_annotated_standardizationAligned.out.sorted.bam
samtools depth -a -m 0 filtered_${name}_annotated_standardizationAligned.out.sorted.bam > filtered_${name}_annotated_standardizationAligned.out.sorted.coverage 
python3 extract_deletions.py filtered_${name}_annotated_standardizationAligned.out.sorted.bam --primer-bed ARTIC_primers_v3.bed --min-deletion-length 20 --ignore-secondary > filtered_${name}_annotated_standardizationAligned.deletion.sorted.txt
```
