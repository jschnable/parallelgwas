#!/usr/bin/bash
#SBATCH --job-name=convertVCFToPLINK
#SBATCH --ntasks=2
#SBATCH --mem-per-cpu=36G
#SBATCH --time=3:05:00
#SBATCH --partition=schnablelab,jclarke,batch,guest
#SBATCH --output=/work/schnablelab/jdavis132/widivGeneExpression/log/%x_%a.out
#SBATCH --error=/work/schnablelab/jdavis132/widivGeneExpression/log/%x_%a.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jdavis132@huskers.unl.edu

ml plink2
ml anaconda
conda activate GWAS

cd $WORK/widivGeneExpression/data

plink2 --vcf twoEnvGenotypes.vcf --make-bed --out twoEnvGenotypes --const-fid
#plink2 --vcf threeEnvGenotypes.vcf --make-bed --out threeEnvGenotypes --const-fid
