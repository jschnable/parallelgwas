#!/usr/bin/bash
#SBATCH --job-name=filterVCF_twoEnvs
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=50G
#SBATCH --time=08:30:00
#SBATCH --partition=schnablelab,jclarke,batch,guest
#SBATCH --output=/work/schnablelab/jdavis132/widivGeneExpression/log/%x_%a.out
#SBATCH --error=/work/schnablelab/jdavis132/widivGeneExpression/log/%x_%a.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jdavis132@huskers.unl.edu

ml vcftools

cd $WORK/widivGeneExpression

vcftools --gzvcf data/WiDiv.vcf.gz --out data/twoEnvGenotypes --keep data/twoEnvGenotypes_wgs.txt --remove-indels --recode
python src/gwas/HomoMAFfiltering.py data/twoEnvGenotypes.recode.vcf data/twoEnvGenotypes.vcf 0.05 0.1
