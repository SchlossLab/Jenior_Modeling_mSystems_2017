Clostridium difficile colonizes alternative nutrient niches during infection across distinct murine gut environments
=======

##Abstract

*Clostridium difficile* infection (CDI) has grown to be the most prevalent cause of hospital acquired infection in the United States. Susceptibility to CDI is induced by recent antibiotic exposure, which is known to alter the structure of the gut microbiome and to affect the availability of growth nutrients in the gut. We hypothesized that *C. difficile* is a generalist that adapts its physiology to the nutrients available within the gut. We orally challenged C57BL/6 mice that previously received one of four antibiotics with *C. difficile* and demonstrated that it was able to colonize the cecum within 18 hours of infection. However, levels of both spore and toxin production, which are known to be affected by nutrient availability, varied between each antibiotic treatment group. To more closely investigate the specific responses of *C. difficile* as it colonized the cecum, we performed *in vivo* transcriptional analysis of *C. difficile* from cecal content of infected mice. This approach revealed variation in expression of genes that drive life-cycle switches as well as metabolic pathways associated with catabolizing a variety of carbon sources such as carbohydrates, amino acids, and amino sugars. To assess which substrates *C. difficile* was most likely exploiting in each antibiotic-perturbed microbiome, we developed a novel metabolite scoring algorithm within the genome-scale bipartite metabolic network of *C. difficile* that incorporated both network topology and transcript abundance to infer the likelihood that a given metabolite was acquired from the environment. Applying this approach, we found that *C. difficile* indeed occupies alternative nutrient niches across each antibiotic-perturbed microbiome and that the highlighted metabolites support significant growth, *in vitro*. Results from this analysis support the hypothesis that consumption of N-acetyl-D-glucosamine and Stickland fermentation substrates are central components of *C. difficile*'s metabolic strategy and pathogenesis. This work has implications for elucidating specifics of the nutrient niche of *C. difficile* during infection and may lead to the discovery of targeted measures to prevent *C. difficile* colonization including potential pre- or probiotic therapies.


Overview
--------
    |- README          # Overview of all content
    |
    |- LICENSE         # Copyright information
    |
    |- Jenior_Modeling_NatMicro_2016.Rmd 	   # executable Rmarkdown for this study
    |
    |- Jenior_Modeling_NatMicro_2016.md 	   # Markdown version of manuscript
    |
    |- manuscript_format.docx 	   # Text tyle-formatting file used for generated docx
    |
    |- nature.csl 	   # Journal-specific formatting csl file
    |
    |- references.bib 	   # Bibtex formatted refereces cited in manuscript
    |
    |- doc/            # Documentation for the study
    |  |- protocols/   # Wet lab and dry lab step-by-step protocols for this study
    |  +- paper/       # Word doc version of manuscript
    |
    |- data            # Raw and primary data
    |  |- README.md         # More specific overview of subdirectory
    |  |- metadata.tsv         # Metadata for mouse experimental groups
    |  |- kegg/  # Reference files used in analysis
    |  |- mapping/         # Normalized read counts mapping to *C. difficile* genes
    |  |- metabolic models/         # Output from bigSMALL analysis
    |  +- wetlab_assays/     # Raw data from wet lab experiments
    |
    |- code/           # Data analysis scripts
    |  |- README.md         # More specific overview of subdirectory
    |  |- R/           # R scipts for figures and analysis
    |      |- figures/      # Code to generatre main body figures
    |      |- supplement/      # Code to generatre supplementary figures
    |      |- support/      # Code to perform additional tasks needed for analysis
    |  |- python/      # python scripts
    |  +- pbs/         # pbs scripts
    |
    |- results         # All output from workflows and analyses
    |  |- README.md         # More specific overview of subdirectory
    |  |- tables/      # Text version of tables (main body)
    |  |- figures/     # Manuscript figures (main body)
    |  +- supplement/  # Supplementary data presentation
    |      |- figures/     # Supplementary figures
    |      +- tables/      # Text and excel versions of supplementary tables
    |
    |- Jenior_Modeling_NatMicro_2016_cache/
    |  +- docx/  # files created during knitting of final docx from Rmd


#### Running the analysis

```
git clone https://github.com/SchlossLab/Jenior_Modeling_NatMicro_2016.git

# Download raw transcriptomic sequencing reads from the SRA (PRJNA354635) to data/fastqs/

qsub code/pbs/db_build.pbs
for transcriptome in cefoperazone_630 clindamycin_630 streptomycin_630 germfree
do
	bash code/pbs/transcriptome.bash transcriptome
done

```
