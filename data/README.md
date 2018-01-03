data/
=======

Overview
--------
    |- README          # Overview of content
    |
    |- metadata.tsv 	   # Metadata for mouse experimental groups
    |
    |- kegg/       # KEGG dictionaries for labeling orthologs and metabolites
    |  |- compound_names.tsv       # KEGG compound code to molecule name dictionary
    |  +- ko_names.tsv     # KO code to enzyme name dictionary
    |
    |- mapping/     # Normalized read counts mapping to *C. difficile* genes
    |  +- cdifficile630
    |     |- all_genes     # reads mapped to all genes in *C. difficile* 630
    |        |- all_infections.ko.RNA.subsample.txt
    |        |- cefoperazone_630.RNA_reads2cdf630.norm.annotated.txt
    |        |- clindamycin_630.RNA_reads2cdf630.norm.annotated.txt
    |        |- streptomycin_630.RNA_reads2cdf630.norm.annotated.txt
    |        +- germfree.RNA_reads2cdf630.norm.annotated.txt
    |     +- select_genes      # mappings to selected genes of interest, not in KEGG
    |        |- cefoperazone_630.RNA_reads2select.all.norm.txt
    |        |- clindamycin_630.RNA_reads2select.all.norm.txt
    |        |- germfree.RNA_reads2select.all.norm.txt
    |        +- streptomycin_630.RNA_reads2select.all.norm.txt
    |
    |- metabolic models/         # Output from bigSMALL analysis
    |  |- combined_sim.tsv     # Shared metabolite importance analysis and confidence intervals
    |  |- optimal_layout.tsv       # Optimum node layout for Fig. 5a
    |  |- expression/      # Raw, and evenly subsamled (.pick) *C. difficil* KO mapping files
    |     |- cefoperazone_630.RNA_reads2cdf630.norm.ko.pick.txt
    |     |- clindamycin_630.RNA_reads2cdf630.norm.ko.pick.txt
    |     |- streptomycin_630.RNA_reads2cdf630.norm.ko.pick.txt
    |     +- germfree.RNA_reads2cdf630.norm.ko.pick.txt
    |  |- cefoperazone_630.bipartite.files/        # bigSMALL output for *C. difficile* during infection of cefoperazone-treated mice
    |     |- bipartite_graph.tsv       # 2 column matrix of all edges in graph
    |     |- cefoperazone_630.importance_score.tsv     # Importance scores and p-values for each metabolite in the nwork
    |     |- cefoperazone_630.mapping.tsv      # Raw KO mapping data from input file
    |     |- cefoperazone_630.topology.tsv     # Topology information for each node in the network
    |     |- compound.lst      # All metabolite nodes
    |     |- enzyme.tst        # All enzyme nodes
    |     |- key_error.log     # Key errors during netowrk annotations
    |     +- parameters.txt        # User defined parameters for importance calculation
    |  |- clindamycin_630.bipartite.files/     # bigSMALL output for *C. difficile* during infection of clindamycin-treated mice (inner files are the same as above)
    |     |- bipartite_graph.tsv
    |     |- clindamycin_630.importance_score.tsv
    |     |- clindamycin_630.mapping.tsv
    |     |- clindamycin_630.topology.tsv
    |     |- compound.lst
    |     |- enzyme.tst
    |     |- key_error.log
    |     +- parameters.txt
    |  |- streptomycin_630.bipartite.files/        # bigSMALL output for *C. difficile* during infection of streptomycin-treated mice (inner files are the same as above)
    |     |- bipartite_graph.tsv
    |     |- streptomycin_630.importance_score.tsv
    |     |- streptomycin_630.mapping.tsv
    |     |- streptomycin_630.topology.tsv
    |     |- compound.lst
    |     |- enzyme.tst
    |     |- key_error.log
    |     +- parameters.txt
    |  +- germfree_630.bipartite.files/        # bigSMALL output for *C. difficile* during infection of germ free mice (inner files are the same as above)
    |     |- bipartite_graph.tsv
    |     |- germfree_630.importance_score.tsv
    |     |- germfree_630.mapping.tsv
    |     |- germfree_630.topology.tsv
    |     |- compound.lst
    |     |- enzyme.tst
    |     |- key_error.log
    |     +- parameters.txt
    |
    |- wetlab_assays/     # Raw data from wet lab experiments
    |  |- cd630_growth.tsv     # Growth curves on important metabolites
    |  |- cef_acetate_630.txt      # GC-MS quantification of acetate in infected mouse ceca
    |  |- cfu.dat      # Colony forming unti quantification from each infection model
    |  +- toxin_titer.dat      # Toxin titer data from each infection model
    |
    |- references/     # Databases for mapping reads to C. difficile genes
    |  |- cdf.kegg.genes.fasta.gz     # Fasta with all C. difficile 630 genes in KEGG (gzipped)
    |  |- cd630.select_genes.fasta.gz      # Fasta with select C. difficile 630 genes from NCBI (gzipped)
    |  |- cdf_db.*.bt2      # Bowtie2 formatted database for cdf.kegg.genes.fasta
    |  +- cd630.select_genes.db.*.bt2      # Bowtie2 formatted database for cd630.select_genes.fasta




