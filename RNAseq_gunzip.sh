#!/bin/bash
#SBATCH --job-name=RNAseq_gunzip 
#SBATCH --time=8:00:00 
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=2 
#SBATCH --mem=2gb 
#SBATCH --mail-type=FAIL
#SBATCH --account=PAS1712

cd $SLURM_SUBMIT_DIR

cd /fs/ess/PAS1712/rcombs/1_pennycress_RNAseq_waterlogging/Workdir/rawdata

gunzip *.fq.gz

#/users/PAS1444/li10917/miniconda2/envs/rnaseq/bin/fastqc -q -t 4 -f fastq *.clean.fq
