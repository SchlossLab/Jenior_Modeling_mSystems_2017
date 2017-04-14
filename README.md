*Clostridium difficile* colonizes alternative nutrient niches during infection across distinct murine gut communities
=======

**Abstract**

*Clostridium difficile* infection (CDI) has become the largest single cause of hospital-acquired infection in the United States. A major risk factor for CDI is prior exposure to antibiotics, as they disrupt the gut bacterial community which protects from *C. difficile* colonization. Multiple antibiotic classes have been associated with susceptibility to CDI, many leading to distinct bacterial communities stemming from variation of bacterial targets. These gut environments present separate metabolic challenges to *C. difficile*, and we hypothesized that *C. difficile* adapts its physiology to available nutrients within differentially susceptible hosts. Utilizing an *in vivo* CDI model, we demonstrated *C. difficile* highly colonized ceca of mice receiving one of three unique antibiotic class pretreatments. We found levels of spores and toxin activity varied between groups, both processes partially regulated by environmental nutrient concentrations. To more closely investigate metabolic responses of *C. difficile* during infection, we performed transcriptomic analysis of the pathogen from the cecal content of infected mice. This revealed hetergeneous expression of numerous catabolic pathways for diverse growth substrates. To assess which resources *C. difficile* was exploiting, we developed a genome-scale metabolic model of *C. difficile* and a transcriptomic-enabled metabolite scoring algorithm that leveraged network architecture. With this platform, we identified nutrients that *C. difficile* used asymmetrically between infection models. These findings were validated through untargeted mass spectrometry of each condition. Our results supported the hypothesis that *C. difficile* inhabits alternative nutrient niches across different gut microbiomes with increased preference for nitrogen-containing carbon sources, specifically Stickland fermentation substrates and host-derived aminoglycans.

**Importance**

Infection by the bacterium *Clostridium difficile* causes an inflammatory diarrheal disease, which can become life-threatening, and has grown to be the most prevalent nosocomial pathogen. Susceptibility to *C. difficile* infection is strongly associated with previous antibiotic treatment, that disrupts the gut microbiota and reduces its ability to prevent colonization. In this study we demonstrate that not only does *C. difficile* alter pathogenesis between hosts pretreated with separate classes of antibiotics, but also exploits different nutrient sources across these environments. Our metabolite importance calculation also provides a platform to study nutrient requirements of pathogens during the context of infection. Our results suggest that *C. difficile* colonization resistance is mediated by multiple groups of bacteria competing for several subsets of nutrients, and could explain why total reintroduction of competitors through fecal microbial transplant is the most effective treatment to date. This work could ultimately contribute to the identification of targeted measures that prevent or reduce *C. difficile* colonization including pre- and probiotic therapies. 



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
