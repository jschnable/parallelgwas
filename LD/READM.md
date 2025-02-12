Steps to calculate LD:

Input vcf filtered file (see the directory to filter vcf data that will help to guide how to filter the vcf data):

Once you have the filtered file vcf data, use this code on terminal to make the bed files or bfiles.

plink2 --vcf example.vcf.gz --make-bed --out ex

Use bfiles calculated to calculate LD for the specific SNP: 

plink --bfile bedfile --r2 --ld-snp SNPname --ld-window-kb 99999 --ld-window 99999 --ld-window-r2 0 --out LD_SNPname


Download the gft file from the phytozome using the gft file generate the gff file 
gffread my.gff3 -T -o my.gtf

Once you have everything:
Run the python script to plot the LD
