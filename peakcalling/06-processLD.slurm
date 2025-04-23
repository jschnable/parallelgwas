#!/usr/bin/bash
#SBATCH --job-name=processLD
#SBATCH --array=1-10
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=3G
#SBATCH --time=00:25:00
#SBATCH --partition=jclarke,schnablelab,batch,guest
#SBATCH --output=/work/schnablelab/jdavis132/widivGeneExpression/log/%x_%a.out
#SBATCH --error=/work/schnablelab/jdavis132/widivGeneExpression/log/%x_%a.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jdavis132@huskers.unl.edu

ml anaconda
conda activate widiv-transcriptome

cd $WORK/widivGeneExpression

R CMD BATCH -data/twoEnvs_chr${SLURM_ARRAY_TASK_ID}.ld src/gwas/processLD.R log/processLD_twoEnvs_${SLURM_ARRAY_TASK_ID}.Rout
R CMD BATCH -data/threeEnvs_chr${SLURM_ARRAY_TASK_ID}.ld src/gwas/processLD.R log/processLD_threeEnvs_${SLURM_ARRAY_TASK_ID}.Rout
