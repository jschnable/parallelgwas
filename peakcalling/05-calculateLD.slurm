#!/usr/bin/bash
#SBATCH --job-name=calculateLD
#SBATCH --array=1-10
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=3G
#SBATCH --time=00:25:00
#SBATCH --partition=jclarke,schnablelab,batch,guest
#SBATCH --output=/work/schnablelab/jdavis132/widivGeneExpression/log/%x_%a.out
#SBATCH --error=/work/schnablelab/jdavis132/widivGeneExpression/log/%x_%a.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jdavis132@huskers.unl.edu

ml plink
ml anaconda
conda activate GWAS

cd $WORK/widivGeneExpression
ld_threshold=0.6

plink --bfile data/twoEnvGenotypes_sigSNPs --r2 --chr ${SLURM_ARRAY_TASK_ID} --ld-window-r2 ${ld_threshold} --ld-window-kb 99999 --out twoEnvs_chr${SLURM_ARRAY_TASK_ID}
plink --bfile data/threeEnvGenotypes_sigSNPs --r2 --chr ${SLURM_ARRAY_TASK_ID} --ld-window-r2 ${ld_threshold} --ld-window-kb 99999 --out threeEnvs_chr${SLURM_ARRAY_TASK_ID}
