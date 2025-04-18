#!/usr/bin/bash
#SBATCH --job-name=calculateEffectiveMarkers
#SBATCH --ntasks=2
#SBATCH --mem-per-cpu=36G
#SBATCH --time=8:05:00
#SBATCH --partition=schnablelab,jclarke,batch,guest
#SBATCH --output=/work/schnablelab/jdavis132/widivGeneExpression/log/%x_%a.out
#SBATCH --error=/work/schnablelab/jdavis132/widivGeneExpression/log/%x_%a.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jdavis132@huskers.unl.edu

ml plink2
ml anaconda
conda activate GWAS

cd $WORK/widivGeneExpression/data

#java -Xmx30g -jar ../src/gec.jar --effect-number --plink-binary twoEnvGenotypes --genome --out twoEnvGenotypes_gec
java -Xmx30g -jar ../src/gec.jar --effect-number --plink-binary threeEnvGenotypes --genome --out threeEnvGenotypes_gec
