#!/bin/bash

## Create Bowtie databases
db_job_id=$(qsub -v transcriptome=$1 code/pbs/db_build.pbs | sed 's/\..*$//')
echo $1 creating mapping databases: $db_job_id

## Assign correct Nextera XT primer pair that was used (5' + 3' reverse complement)
F_primer=AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT
if [ $1 = 'cefoperazone_630' ]; then
	R_primer=CAAGCAGAAGACGGCATACGAGATCGTGATGTGACTGGAGTTCAGACGTGTGCTCTTCCGATC
	infection='infected'
elif [ $1 = 'clindamycin_630' ]; then
	R_primer=CAAGCAGAAGACGGCATACGAGATTGGTCAGTGACTGGAGTTCAGACGTGTGCTCTTCCGATC
	infection='infected'
elif [ $1 = 'streptomycin_630' ]; then
	R_primer=CAAGCAGAAGACGGCATACGAGATAAGCTAGTGACTGGAGTTCAGACGTGTGCTCTTCCGATC
	infection='infected'
elif [ $1 = 'germfree' ]; then
	R_primer=CAAGCAGAAGACGGCATACGAGATATCCACTCGTGACTGGAGTTCAGACGTGTGCTCTTCCGATC
	infection='infected'
fi

## Pool reads from the same lane
pool_job_id=$(qsub -v transcriptome=$1 code/pbs/pool.pbs | sed 's/\..*$//')
echo $1 read pooling: $pool_job_id

## Curate reads
trimming_job_id=$(qsub -v transcriptome=$1,forward=$F_primer,reverse=$R_primer -W depend=afterok:$pool_job_id code/pbs/trimming.pbs | sed 's/\..*$//')
echo $1 quality trimming: $trimming_job_id

## Map curated reads to C. difficile
mapping_job1_id=$(qsub -v transcriptome=$1 -W depend=afterok:$filter_job_id code/pbs/map2select_630.pbs | sed 's/\..*$//')
echo $1 mapping job 1 reads: $mapping_job1_id
mapping_job2_id=$(qsub -v transcriptome=$1 -W depend=afterok:$filter_job_id code/pbs/map2cdf_genes.pbs | sed 's/\..*$//')
echo $1 mapping job 1 reads: $mapping_job2_id


echo