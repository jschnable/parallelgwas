# parallelgwas
Code to help in running and analyzing eQTL and FarmCPU GWAS results across large numbers of traits.

## Suggested workflow for parallel GWAS: 
1. Create conda environment GWAS with R packages tidyverse and rMVP installed.
2. Filter VCF file using ``00-filterVCF.slurm``. Requires list of genotypes to keep as a text file.
3. Convert filtered VCF to PLINK and rMVP formats using ``01-xxx`` SLURM scripts.
4. Estimate effective SNP number using GEC (``02-estimateEffectiveSNPs.slurm``).
5. Split phenotype files using ``splitDataFrame.R`` or create sample files for resampling GWAS using ``samplePhenotypesForResampling.R``.
6. Create model specification text files:
   - pheno_file_list: path to each phenotype file from step 5
   - genostem_list: if using multiple different genotype files, path to genotype file each phenotype should be run with, in the same order. If using one genotype file, this path can be specified as a static value in `g` in ``03-runXX.slurm``
   - totalmarker_list: See genostem notes. Should contain total number of total markers in corresponding genotype file.
   - effectivemarker_list: See genostem notes. Should contain effective number of SNPs (generated in output from step 4) in corresponding genotype file.
7. Run GWAS! Change paths in 03-runXX to match files created in step 6.
8. Summarise output files using ``summariseSignals.R`` (MLM) or ``getSupport.R`` (FarmCPU).

## Suggested workflow for peak calling based on LD. All scripts are in peakcalling/
1. Generate list of all SNPs above significance threshold in at least one GWAS model. Save in a text file. 
2. Filter PLINK file used as input for GEC to significant SNPs only using ``04-filterPLINKToSigSNPs.slurm``.
3. Calculate LD within a chromosome between significant SNPs using ``05-calculateLD.slurm``.
4. Process plink-generated LD files for input to peak calling script using ``06-processLD.slurm``.
5. Generate text files listing the  traits to process in each job of the array. Each list should only have traits that will use the same LD file to call peaks (i.e., all traits have at least the minimum number of SNPs for a peak on chromosome 1, etc). Example code for generating these files is in ``generatePeakCallingTextFiles.R``. If calling peaks to collapse FarmCPU resampling signals across environments for the same trait, set ``min_snps_per_peak=1``. Also generate the following text files (maintaining order of lines across files): 
   - ld_file_list: path to ld file to use 
   - signal_file_list: path to signal file to use
   - traitList_file_list: path to trait list file (see above) to use
   - chrom_list: chromosome being processed 
6. Call peaks! Change paths in ``07-callPeaks.slurm`` to match files created in step 5 and run. 
7. Summarise output files using ``summarisePeaks.R``.

