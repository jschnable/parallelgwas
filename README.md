# parallelgwas
Code to help in running and analyzing eQTL and FarmCPU GWAS results across large numbers of traits.

# Suggested workflow: 
1. Create conda environment GWAS with R packages tidyverse and rMVP installed.
2. Filter VCF file using 00-filterVCF.slurm. Requires list of genotypes to keep as a text file.
3. Convert filtered VCF to PLINK and rMVP formats using 01-xxx SLURM scripts.
4. Estimate effective SNP number using GEC (02-estimateEffectiveSNPs.slurm).
5. Split phenotype files using splitDataFrame.R or create sample files for resampling GWAS using samplePhenotypesForResampling.R.
6. Create model specification text files:
   - pheno_file_list: path to each phenotype file from step 5
   - genostem_list: if using multiple different genotype files, path to genotype file each phenotype should be run with, in the same order. If using one genotype file, this path can be specified as a static value in g in 03-runXX.slurm
   - totalmarker_list: See genostem notes. Should contain total number of total markers in corresponding genotype file.
   - effectivemarker_list: See genostem notes. Should contain effective number of SNPs (generated in output from step 4) in corresponding genotype file.
  

  7. Run GWAS! Change paths in 03-runXX to match files created in step 6.
