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
