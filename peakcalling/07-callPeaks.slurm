#!/usr/bin/bash
#SBATCH --array=1-306
#SBATCH --job-name=callPeaks
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=15
#SBATCH --mem-per-cpu=5G
#SBATCH --time=6:00:00
#SBATCH --partition=schnablelab,jclarke,batch,guest
#SBATCH --output=/work/schnablelab/jdavis132/widivGeneExpression/log/%x_%a.out
#SBATCH --error=/work/schnablelab/jdavis132/widivGeneExpression/log/%x_%a.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jdavis132@huskers.unl.edu

cd /work/schnablelab/jdavis132/widivGeneExpression
ml anaconda
conda activate widiv-transcriptome

ld_file_list="data/peakCalling_ld_file_list.txt"
signal_file_list="data/peakCalling_signal_file_list.txt"
traitList_file_list="data/peakCalling_traitList_file_list.txt"
chrom_list="data/peakCalling_chrom_list.txt"

ld=`sed -n "${SLURM_ARRAY_TASK_ID}"p ${ld_file_list} |  awk '{print $1}'`
signal=`sed -n "${SLURM_ARRAY_TASK_ID}"p ${signal_file_list} |  awk '{print $1}'`
traitListFile=`sed -n "${SLURM_ARRAY_TASK_ID}"p ${traitList_file_list} |  awk '{print $1}'`
chrom=`sed -n "${SLURM_ARRAY_TASK_ID}"p ${chrom_list} |  awk '{print $1}'`

R CMD BATCH -data/${ld} -data/${signal} -${chrom} -data/${traitListFile} -${SLURM_ARRAY_TASK_ID} src/gwas/callPeaks.R log/callPeaks_${SLURM_ARRAY_TASK_ID}.Rout
