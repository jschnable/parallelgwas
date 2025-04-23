#!/usr/bin/bash
#SBATCH --job-name=filterPLINK
#SBATCH --ntasks=2
#SBATCH --mem-per-cpu=30G
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

plink --bfile data/twoEnvGenotypes --extract data/sigSNPs_twoEnvs.txt --make-bed --out data/twoEnvGenotypes_sigSNPs
#plink --bfile data/threeEnvGenotypes --extract data/sigSNPs_threeEnvs.txt --make-bed --out data/threeEnvGenotypes_sigSNPs
