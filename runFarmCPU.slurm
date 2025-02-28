#!/usr/bin/bash
#SBATCH --array=1-80
#SBATCH --job-name=farmcpu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5
#SBATCH --mem-per-cpu=15G
#SBATCH --time=01:45:00
#SBATCH --partition=schnablelab,jclarke,batch,guest
#SBATCH --output=/work/schnablelab/jdavis132/widivGeneExpression/log/%x_%a.out
#SBATCH --error=/work/schnablelab/jdavis132/widivGeneExpression/log/%x_%a.err

#change this path
cd $WORK/widivGeneExpression
ml anaconda
conda activate GWAS # conda environment with rMVP installed

# ensure lines in these files are in the same order as each other
pheno_file_list="data/gwas_phenofiles.txt" # file with path to phenotype files to run
genostem_list="data/gwas_genostems.txt" # file with filestems of rMVP-formatted genotype data
totalmarker_list="data/gwas_totalmarkers.txt" # file with number of total markers in the respective genotype files
effectivemarker_list="data/gwas_effectivemarkers.txt" # file with number of effective markers in respective genotype file

p=`sed -n "${SLURM_ARRAY_TASK_ID}"p ${pheno_file_list} |  awk '{print $1}'`
g=`sed -n "${SLURM_ARRAY_TASK_ID}"p ${genostem_list} |  awk '{print $1}'`
total_markers=`sed -n "${SLURM_ARRAY_TASK_ID}"p ${totalmarker_list} |  awk '{print $1}'`
effective_markers=`sed -n "${SLURM_ARRAY_TASK_ID}"p ${effectivemarker_list} |  awk '{print $1}'`

cd $SCRATCH
# change these paths
cp /work/schnablelab/jdavis132/widivGeneExpression/${p} .
cp /work/schnablelab/jdavis132/widivGeneExpression/${g}.* .
cp /work/schnablelab/jdavis132/widivGeneExpression/src/gwas/FarmCPU.R .

R CMD BATCH -${p} -${g} -${total_markers} -${effective_markers} FarmCPU.R FarmCPU_${SLURM_ARRAY_TASK_ID}.Rout

# change these paths
cp *signals.csv /work/schnablelab/jdavis132/widivGeneExpression/out/gwas/${g}/
cp *.Rout /work/schnablelab/jdavis132/widivGeneExpression/log/
