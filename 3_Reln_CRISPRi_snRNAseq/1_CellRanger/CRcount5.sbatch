#!/bin/bash

#SBATCH -J CRcount5_2023-JD-0062
#SBATCH -o %x.out
#SBATCH -e %x.err
#SBATCH -t 48:00:00
#SBATCH -p medium
#SBATCH -c 12
#SBATCH --mem-per-cpu=12G
#SBATCH --mail-user=cenewman@uab.edu
#SBATCH --mail-type=FAIL,END

export DATA_DIR=/data/project/daylab/2023-JD-0062/
export WORK_DIR=/scratch/cenewman/CRcount5
cd $WORK_DIR

module load CellRanger

cellranger count --id=CR_outputs_2023-JD-0062 \
                 --fastqs=$DATA_DIR/fastqs/Reln_S2 \
                 --transcriptome=/data/user/cenewman/reference/mRatBN7_2 \
                 --localcores=12 \
                 --localmem=144

date

exit 0
