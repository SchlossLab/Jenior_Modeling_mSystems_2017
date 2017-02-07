*Clostridium difficile* colonizes alternative nutrient niches during infection across distinct murine gut environments
=======

###Abstract

*Clostridium difficile* infection (CDI) has become the largest single cause of hospital-acquired infection in the United States. A compromised gut microbiota, typically through recent antibiotic exposure, is a prerequisite feature of *C. difficile* colonization susceptibility. This has been described for multiple antibiotic classes in which many result in distinct gut communities, each presenting individual metabolic challenges to *C. difficile*. We hypothesized that *C. difficile* must adapt its physiology to nutrients availability within the gut. Utilizing an *in vivo* model of CDI, we demonstrated that *C. difficile* highly colonized the cecum of mice that received one of three antibiotic pretretments. We found levels of spore and toxin production varied between each antibiotic treatment group, and both processes are knwon to be regulated by specific nutrient concnetrations. To more closely investigate specific responses of *C. difficile* during infection, we performed transcriptional analysis of *C. difficile* from cecal content of infected mice. This revealed variation in expression of life-cycle switches and catabolic pathways for a variety of carbon sources. In order to assess which substrates *C. difficile* was exploiting, we further characterized the systems with transcriptomic-enabled genome-scale metabolic modeling and untargeted metabolomic analysis. Through the development of a novel metabolite scoring algorithm, leveraging the metabolic model architecture, we were able to infer that a given metabolite was acquired from the environment. Our results support the hypothesis that *C. difficile* indeed occupies alternative nutrient niches by metabolizing separate carbohydrate sources in each infection and these distinctions track with disparity seen in pathogenicity. Additionally, these data highlight conserved elements of *C. difficile*'s metabolic strategy across infections, including the consumption of N-acetyl-D-glucosamine and Stickland fermentation substrates. 

###Importance

In this study we demonstrate that not only does *C. difficile* alter pathogenesis between differentially sensitized hosts, but also exploits separate nutrient niches across environments. Our results support that *C. difficile* possesses a highly plastic nutrient niche space, allowing it to successfully infect distinct hosts and ultimately cause disease. This work also provides evidence that *C. difficile* virulence may be driven by accessibility of specific carbohydrates utilized for growth during each infection. This work has implications for elucidating drivers of *C. difficile* pathogenesis and uncover specifics colonization resistance. This could lead to the discovery of targeted measures to prevent *C. difficile* colonization including potential pre- or probiotic therapies. Furthermore, the metabolite importance calculation workflow descibed here could provide a useful platform to enable more rapid discoveries for the nutrient requirements of bacteria to be made in the future.


Overview
--------
    |- README          # Overview of all content
    |
    |- LICENSE         # Copyright information
    |
    |- Jenior_Modeling_mBio_2016.Rmd 	   # executable Rmarkdown for this study
    |
    |- Jenior_Modeling_mBio_2016.md 	   # Markdown version of manuscript
    |
    |- manuscript_format.docx 	   # Text tyle-formatting file used for generated docx
    |
    |- mbio.csl 	   # Journal-specific formatting csl file
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
    |- Jenior_Modeling_mBio_2016_cache/
    |  +- docx/  # files created during knitting of final docx from Rmd


#### Running the analysis

```
# Clone bigSMALL to Desktop
git clone git@github.com:mjenior/bigsmall.git
git clone git@github.com:SchlossLab/Jenior_Modeling_NatMicro_2016.git

# Download raw transcriptomic sequencing reads from the SRA (PRJNA354635) to data/fastqs/

qsub code/pbs/db_build.pbs
for transcriptome in cefoperazone_630 clindamycin_630 streptomycin_630 germfree
do
	bash code/pbs/transcriptome.bash transcriptome
done

python ~/Desktop/bigsmall/bigsmall.py data/metabolic_models/expression/streptomycin_630.RNA_reads2cdf630.norm.ko.pick.txt --iters 10000 --name streptomycin_630
python ~/Desktop/bigsmall/bigsmall.py data/metabolic_models/expression/cefoperazone_630.RNA_reads2cdf630.norm.ko.pick.txt --iters 10000 --name cefoperazone_630
python ~/Desktop/bigsmall/bigsmall.py data/metabolic_models/expression/clindamycin_630.RNA_reads2cdf630.norm.ko.pick.txt --iters 10000 --name clindamycin_630
python ~/Desktop/bigsmall/bigsmall.py data/metabolic_models/expression/germfree.RNA_reads2cdf630.norm.ko.pick.txt --iters 10000 --name germfree_630

R code/R/figures/fig_1.R
R code/R/figures/fig_2.R
R code/R/figures/fig_3.R
R code/R/figures/fig_4.R
R code/R/figures/fig_5.R
R code/R/figures/fig_6.R
R code/R/figures/fig_7.R

R code/R/supplement/fig_S1.R
R code/R/supplement/fig_S2.R
R code/R/supplement/fig_S3.R
R code/R/supplement/fig_S4.R
R code/R/supplement/fig_S5.R

```
