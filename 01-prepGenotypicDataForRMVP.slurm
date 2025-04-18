#!/usr/bin/bash
#SBATCH --job-name=prepGenotypicData
#SBATCH --ntasks=2
#SBATCH --mem-per-cpu=36G
#SBATCH --time=12:00:00
#SBATCH --partition=schnablelab,jclarke,batch,guest
#SBATCH --output=/work/schnablelab/jdavis132/widivGeneExpression/log/%x_%a.out
#SBATCH --error=/work/schnablelab/jdavis132/widivGeneExpression/log/%x_%a.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jdavis132@huskers.unl.edu

ml anaconda
conda activate GWAS

cd $WORK/widivGeneExpression

#R CMD BATCH -data/twoEnvGenotypes.vcf -data/twoEnvGenotypes.mvp src/gwas/prepGenotypicDataForRMVP.R
R CMD BATCH -data/threeEnvGenotypes.vcf -data/threeEnvGenotypes.mvp src/gwas/prepGenotypicDataForRMVP.R
