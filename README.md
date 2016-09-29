Colonization Resistance Metatranscriptomics
=======

Research project initialization and organization following reproducible research
guidelines as modified for use with typical microbial ecology projects

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
