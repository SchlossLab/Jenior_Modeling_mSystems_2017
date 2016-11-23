code/
=======

Overview
--------
    |- README          # Overview of content
    |
    |- python/
    |  |- read_growth.py       # Formats raw growth meter data for analysis
    |  |- count_idxstats.py       # Counts the number of reads the mapped in an idxstats file
    |  |- normalize_idxstats.py       # Normalizes mapping abundances to gene and read length
    |  +- seq_stats.py     # Reports various metrics on input fasta/fastq files
    |
    |- R/     # R scripts for analysis and figure generation
    |  |- figures/     # Code to generate main body figures
    |     |- fig_1.R
    |     |- fig_2.R
    |     |- fig_3.R
    |     |- fig_4.R
    |     |- fig_5.R
    |     +- fig_6.R
    |  |- supplement/     # Code for supplementary figures
    |     |- fig_S1.R
    |     |- fig_S2.R
    |     +- fig_S3.R
    |  +- support/     # Additional analysis scripts
    |     |- curve_test.R     # Formates growth data for ANOVA tests
    |     +- subsample.R     # Evenly subsamples mapping files for even analysis
    |
    |- pbs/     # PBS scripts used during data processing
    |  |- db_build.pbs     # Creates bowtie databases
    |  |- pool.pbs      # Pools raw reads
    |  |- trimming.pbs      # Quality and adapter trims pooled reads
    |  |- map2cdf_genes.pbs      # Map to all C. difficile 630 genes
    |  |- map2select_630.pbs      # Map to specific C. difficile 630 genes of interest
    |  +- transcriptome.bash      # Bash script to exectue whole pipeline
