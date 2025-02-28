#!/usr/bin/bash
#SBATCH --job-name=prepGenotypicData
#SBATCH --ntasks=2
#SBATCH --mem-per-cpu=36G
#SBATCH --time=12:00:00
#SBATCH --partition=schnablelab,jclarke,batch,guest
#SBATCH --output=/work/schnablelab/jdavis132/widivGeneExpression/log/%x_%a.out
#SBATCH --error=/work/schnablelab/jdavis132/widivGeneExpression/log/%x_%a.err

ml anaconda
conda activate GWAS # conda environment with rMVP installed
# change this path
cd $WORK/widivGeneExpression
# change these paths; first arg is path to vcf file to prep, second is filestem for output genotype files, third is path to script, 4th is path to log file to write R console output
R CMD BATCH -data/twoEnvGenotypes.vcf -twoEnvGenotypes.mvp src/gwas/prepGenotypicDataForRMVP.R
R CMD BATCH -data/threeEnvGenotypes.vcf -threeEnvGenotypes.mvp src/gwas/prepGenotypicDataForRMVP.R
