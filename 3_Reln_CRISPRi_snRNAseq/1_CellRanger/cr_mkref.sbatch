#!/bin/bash

#SBATCH -J CR_mkref
#SBATCH -o %x.out
#SBATCH -e %x.err
#SBATCH -t 12:00:00
#SBATCH -p short
#SBATCH -c 6
#SBATCH --mem-per-cpu=16G
#SBATCH --mail-user=cenewman@uab.edu
#SBATCH --mail-type=FAIL,END

export DATA_DIR=/data/project/daylab/2023-JD-0062
export WORK_DIR=/scratch/cenewman/mkref
cd $WORK_DIR

module load CellRanger

cellranger mkref --genome=mRatBN7_2 \
                 --fasta=$DATA_DIR/custom_rn7_outs/Rattus_norvegicus.mRatBN7.2.dna.toplevel_Chimera.fa \
                 --genes=$DATA_DIR/custom_rn7_outs/GTF/Rattus_norvegicus.mRatBN7.2.109_Chimera.gtf
                 --nthreads=6

date

exit 0
