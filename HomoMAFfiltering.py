from collections import Counter
import sys


pathtovcf=sys.argv[1]
pathtofilteredvcf=sys.argv[2]
mafthrehsold=sys.argv[3]
hetthreshold=sys.argv[4]

data=open(str(pathtovcf))

mafthreshold=float(mafthrehsold)
hetthreshold=float(hetthreshold)
with open(str(pathtofilteredvcf), 'w') as SNPs:
    for line in data:
        line=line.strip()
        if line.startswith("#"):
            # print(line)
            SNPs.write(f'{line}\n')
            continue

        fields=line.split("\t")
        # print(fields)

        pos=fields[2]

        genotypes=fields[9:]

        genotype_counts = Counter(genotypes)
        hom_ref = genotype_counts.get('0/0',0) + genotype_counts.get('0|0',0)
        hom_alt = genotype_counts.get('1/1',0) + genotype_counts.get('1|1',0)
        het_count = genotype_counts.get('0/1',0) + genotype_counts.get('0|1',0) +genotype_counts.get('1/0',0) + genotype_counts.get('1|0',0)
        total_count=len(genotypes)

        hom_maf = min(hom_ref, hom_alt) / total_count
        het_f = het_count/total_count

        if hom_maf > mafthreshold and het_f < hetthreshold:
                SNPs.write(f'{line}\n')
