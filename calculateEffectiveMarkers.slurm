#!/usr/bin/bash
#SBATCH --job-name=calculateEffectiveMarkers
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=6G
#SBATCH --time=2:05:00
#SBATCH --partition=jclarke,batch,guest
#SBATCH --output=/work/schnablelab/user/project/log/%x_%a.out
#SBATCH --error=/work/schnablelab/user/project/log/%x_%a.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=user@huskers.unl.edu

# change the user, project dir, and vcf_stem values above & in next 2 lines
PROJECT_DIR=/work/schnablelab/user/project
VCF_STEM=project
ml plink2
cd ${PROJECT_DIR}
# convert vcf to plink
plink2 --vcf ${VCF_STEM}.vcf --make-bed --out ${VCF_STEM}
# estimate effective SNP number 
java -Xmx5g -jar gec.jar --effect-number --plink-binary ${VCF_STEM} --genome --out ${VCF_STEM}_gec
# results in ${PROJECT_DIR}/${VCF_STEM}_gec.sum
