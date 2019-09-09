#!/bin/bash -e

# Job Name
#$ -N r3_d28

# Execute the script from the Current Working Directory
#$ -cwd

# Merge the output of the script, and any error messages generated to one file
#$ -j y

# Send the output of the script to a directory called 'UGE-output' uder current
# working directory (cwd)
if [ ! -d "UGE-output" ]; then #Create output directory in case it does NOT
          exist
mkdir UGE-output
fi
#$ -o UGE-output/

# Tell the job your hardware requirements
#$ -l mem_free=64G

# Send mail when the job is submitted, and when the job
# completes
#$ -m be

#  Specify an email address to use
#$ -M rohit.farmer@.nih.gov

source ~/condarc
conda activate h5n1
Rscript --vanilla SCRIPTS/eNetXplorer/eNetXplorer_R3_180530_d28.R


