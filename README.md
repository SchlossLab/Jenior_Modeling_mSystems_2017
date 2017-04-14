*Clostridium difficile* colonizes alternative nutrient niches during infection across distinct murine gut communities
=======

**Abstract**

*Clostridium difficile* infection (CDI) has become the largest single cause of hospital-acquired infection in the United States. CDI susceptibility is most frequently associated with previous antibiotic exposure, which disrupt the gut bacterial community. This has been described for multiple antibiotic classes which result in distinct gut communities, each presenting separate metabolic challenges to *C. difficile*. We hypothesized that *C. difficile* adapts its physiology to nutrient availability within differentially susceptible gut environments. Utilizing an *in vivo* model of CDI, we demonstrated that *C. difficile* highly colonized the cecum of mice receiving one of three individual antibiotic pretreatments. We found levels of spore and toxin production varied between each group, both processes partially regulated by environmental nutrient concentrations. To more closely investigate metabolic responses of *C. difficile* during infection, we performed transcriptomic analysis of the pathogen from cecal content of infected mice. This revealed expression variation in numerous catabolic pathways for various carbon sources. To assess which substrates *C. difficile* was exploiting, we developed a transcriptomic-enabled genome-scale metabolic model of *C. difficile* and a metabolite scoring algorithm that leveraged network architecture. With this platform, we identified carbon sources used by *C. difficile* asymmetrically between infection models. These results were validated through correlation with untargeted mass spectrometry analysis from each condition. Our results supported the hypothesis that *C. difficile* indeed metabolized alternative carbon sources across colonized environments. These data also highlighted conserved elements of *C. difficile*'s metabolic strategy, specifically consumption of host-derived aminoglycans and Stickland fermentation substrates.

**Importance**

In this study we demonstrate that not only does *C. difficile* alter pathogenesis between differentially sensitized hosts, but also exploits separate nutrient niches across environments. Our results support that *C. difficile* possesses a highly plastic nutrient niche space, allowing it to successfully infect distinct communities. This work also provides evidence that *C. difficile* virulence may be driven by accessibility of specific carbon sources during infection. This work could lead to the discovery of targeted measures to prevent *C. difficile* colonization including potential pre- or probiotic therapies. Our metabolite importance calculation workflow also provides a platform to the study of nutrient requirements of pathogens in the context of infection or even patterns of substrate utilization in communities of bacteria.



Overview
--------
    |- README          # Overview of all content
    |
    |- LICENSE         # Copyright information
    |
    |- Jenior_Modeling_mBio_2017.Rmd 	   # executable Rmarkdown for this study
    |
    |- Jenior_Modeling_mBio_2017.md 	   # Markdown version of manuscript
    |
    |- manuscript_format.docx 	   # Text tyle-formatting file used for generated docx
    |
    |- mbio.csl 	   # Journal-specific formatting csl file
    |
    |- references.bib 	   # Bibtex formatted refereces cited in manuscript
    |
    |- protocols/            # Wet lab and dry lab step-by-step protocols for this study
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
    |      +- tables/      # Excel versions of supplementary tables


#### Running the analysis

```
# Clone bigSMALL to Desktop
git clone git@github.com:mjenior/bigsmall.git
git clone git@github.com:SchlossLab/Jenior_Modeling_MBio_2017.git

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
R code/R/supplement/fig_S1.R
R code/R/supplement/fig_S2.R
R code/R/supplement/fig_S3.R
R code/R/supplement/fig_S4.R
R code/R/supplement/fig_S5.R
R code/R/supplement/fig_S6.R

```
