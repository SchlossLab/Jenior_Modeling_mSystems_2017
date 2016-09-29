Clostridium difficile colonizes alternative nutrient niches during infection across unique murine gut environments
=======

Abstract
--------

	Clostridium difficile infection (CDI) as grown to be the greatest cause of hospital acquired infection in the United States. Susceptibility to CDI is induced by previous antibiotic exposure, which has been shown to alter the structure of the gut microbiome. These changes have been associated with changes in bacterial growth nutrient availability in the gut, often increasing concentrations of several useable by C. difficile. In this study, we orally challenged C57BL/6 mice with C. difficile str. 630 and demonstrated that it was able to colonize the ceca in three separate models of antibiotic induced susceptibility to the same high degree (~1×108 CFU/g content) within 18 hours of inoculation. However, despite equal vegetative cell load at this time point, the levels of both spore and toxin production vary between each antibiotic treatment group. The expression of both phenotypes have both been linked to environmental concentrations of certain substrates, and this indicated possible differences in the nutrient niche C. difficile inhabits between susceptible gut conditions. To more closely investigate the specific responses of C. difficile as it colonizes the GI tract of mice, we performed in vivo C. difficile-focused RNA-Seq analysis from cecal content of infected mice. This approach identified differences in expression for genes associated with life-cycle stages and pathogenicity between antibiotic pretreatments, in agreement with previous results. We then went on to observe numerous variations between condition in metabolic pathways associated with carbohydrate and amino acid catabolism, indicating that C. difficile likely colonizes alternative nutrient niches across the environments it colonizes. In order to assess which aspects of the gut environment C. difficile is exploiting during infection, we sought to identify the growth nutrients that are most likely being used the the pathogen in each colonized condition. To accomplish this we developed a novel substrate scoring algorithm within the genome-scale bipartite metabolic network of C. difficile...

Overview
--------

    project
    |- README          # the top level description of content
    |
    |- doc/            # documentation for the study
    |  |- protocols/   # wetlab and drylab
    |  +- paper/       # manuscript
    |
    |- data            # raw and primary data
    |  |- references/  # reference files  used in analysis
    |  |- raw/         # raw data, unaltered
    |  +- process/     # cleaned data
    |
    |- code/           # data analysis scripts
    |  |- R/           # R
    |  |- python/      # python
    |  +- pbs/         # pbs
    |
    |- results         # all output from workflows and analyses
    |  |- tables/      # text version of tables to be rendered with kable in R
    |  |- figures/     # manuscript figures
    |  +- supplement/  # supplementary data presentation
    |  |   |- figures/     # figures
    |  |   +- tables/      # text version of tables
    |
    |- study.Rmd       # executable Rmarkdown for this study, if applicable
    |- study.md        # Markdown (GitHub) version of the *Rmd file
