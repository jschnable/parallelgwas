import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
import seaborn as sns
import numpy as np
import re
import pandas as pd
import matplotlib
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
from matplotlib import cm, colors
from matplotlib.colors import ListedColormap

#configurable variables
target_chr = 'Chr09'
ld_file = "output.ld" #(write the name of the .ld file generated)
target_loc = 62620720 #(write the position of the sap)
windowsize1_left = 30000
windowsize2_right = 30000
stepsize = 1000
myfigsize = (3,1.5)
xlab_string = "Position on Chromosome 9" #(fix this as per required)
savename = "Test_LD_plot"
candidate_gene = "Sobic.009G254200" #(name of the gene)


plt.rcParams.update({
    "text.usetex": True,
    "font.family": "sans-serif",
    "font.sans-serif": "Helvetica",
    "font.weight": "normal",
    "font.size": 6
})

def ld_parse(afile):
    fh = open(afile)
    fh.readline()
    ld_dict = {}
    for x in fh:
        y = x.split()
        # print(y)
        ld_dict[y[-2]] = [float(y[-1]), y[-3]]
    return ld_dict

def genomicpos(x,pos):
    x = x/1000000
    return "{:,.2f}".format(x)

def gff2strand(afile):
    fh = open(afile)
    stranddict = {}
    for x in fh:
        if x[0] == '#': continue
        y = x.strip().split('\t')
        if y[2] != 'gene': continue
        try:
            mychr = str(int(y[0].replace('Chr','')))
        except:
            continue
        mystrand = y[6]
        mystart = int(y[3])
        mystop = int(y[4])
        myname = myname = y[-1].split('.v3')[0].replace('ID=','')
        stranddict[myname] = mystrand
    return stranddict


def defparse(astr):
    dd = {}
    y = astr.strip().split(';')
    for a in y:
        b = a.split('=')
        dd[b[0]] = b[1]
    return dd

def gff_parse(afile,tchr,rstart,rstop):
    fh = open(afile)
    gene_dict = {}
    for x in fh:
        if x[0] == '#': continue
        y = x.strip().split('\t')
        if y[0] != "{0}".format(str(tchr).zfill(1)): continue
        if int(y[3]) > rstop: continue
        if int(y[4]) < rstart: continue
        defdict = defparse(y[-1])
        if y[2] == 'mRNA' and defdict["longest"] and defdict["longest"] == '1':
            mygene = defdict["ID"]
            gene_dict[mygene] = {"start":int(y[3]),"stop":int(y[4]),"utr":[],"cds":[],"strand":y[6]}
        if y[2] == 'CDS':
            mygene = defdict["Parent"]
            if not mygene in gene_dict: continue
            gene_dict[mygene]["cds"].append((int(y[3]),int(y[4])-int(y[3])))
        if y[2] == "five_prime_UTR" or y[2] == "five_prime_UTR":
            mygene = defdict["Parent"]
            if not mygene in gene_dict: continue
            gene_dict[mygene]["utr"].append((int(y[3]),int(y[4])-int(y[3])))
    return(gene_dict)
    
pter = FuncFormatter(genomicpos)

colors = ['#808080','#FF0000', '#0000FF']  # Red, Blue, Grey
cmap = ListedColormap(colors)

gene_dict = gff_parse("Sbicolor_730_v5.1.gene_exons.gff3",target_chr,target_loc-windowsize1_left,target_loc+windowsize2_right) #name of the .gff3 file change as per needed

ld_dict = ld_parse(ld_file)

df=pd.DataFrame.from_dict(ld_dict).T
df.columns=['r2', 'position']
df=df.sort_values('position', ascending=False)
df['SNP']=df.index
df=df.reset_index()

df.position=pd.to_numeric(df.position)
df.r2=pd.to_numeric(df.r2)
df = df[(df['position'] >= target_loc - windowsize1_left) & (df['position'] <= target_loc + windowsize2_right)]

fig=plt.figure(figsize=myfigsize)
ax=fig.add_subplot(1,1,1)
#norm = plt.Normalize(df['r2'].min(), df['r2'].max())
#sm = plt.cm.ScalarMappable(cmap="flare", norm=norm)
#sm.set_array([])

ax=sns.scatterplot(data=df, x='position', y='r2', hue='r2',hue_norm=(0,1), palette="flare", legend=False,s=50)

# Add colorbar
cbar = fig.colorbar(cm.ScalarMappable(cmap=matplotlib.colormaps["flare"]), ax=ax)

#cbar=plt.colorbar(sm, ax=ax)
#sbarticks=cbar.get_ticks()
#cbar.set_ticks(sbarticks)
#cbar.set_ticklabels([round(c,2) for c in sbarticks])
#cbar.set_label('$R^2$',labelpad=-70, y=0.5)

plt.ylim(-0.25,1.1)


for g in gene_dict:
    # print(g)
    # myname = g.split('_')[0].split('\.1|\.2|\.3|\.4|\.5')[0]
    myname=re.split(r'\.1|\.2|\.3|\.4|\.5', g)[0]
    # ax.broken_barh(gene_dict[g]["utr"],(-0.3,.1),color="blue")
    # ax.broken_barh(gene_dict[g]["cds"],(-0.3,.1),color="green")
    if gene_dict[g]["strand"] == '+':
        arrow = mpatches.Arrow(gene_dict[g]["start"],-0.15,gene_dict[g]["stop"]-gene_dict[g]["start"]+100,0,width=0.1,color="black")
    if gene_dict[g]["strand"] == '-':
        arrow = mpatches.Arrow(gene_dict[g]["stop"],-0.15,gene_dict[g]["start"]-gene_dict[g]["stop"]+100,0,width=0.1,color="black")
    ax.broken_barh(gene_dict[g]["utr"],(-0.2,.1),color="blue")
    ax.broken_barh(gene_dict[g]["cds"],(-0.2,.1),color="green")
    ax.add_patch(arrow)
    # print(myname)
    # if (gene_dict[g]["start"]+gene_dict[g]["stop"])/2 < target_loc - windowsize or (gene_dict[g]["start"]+gene_dict[g]["stop"])/2  > target_loc + windowsize: continue
    if myname==candidate_gene:
        # print(myname)
        ax.text((gene_dict[g]["start"]+gene_dict[g]["stop"])/2,-0.05,myname,ha="center")
        ax.text((gene_dict[g]["start"]+gene_dict[g]["stop"])/2,-0.08,myname+' (Y1)',ha="center")
        ax.broken_barh(gene_dict[g]["utr"],(-0.2,.1),color="blue")
        ax.broken_barh(gene_dict[g]["cds"],(-0.2,.1),color="red")
#yticks=ax.get_yticks()
#ax.set_yticks(yticks[2:8])
#ax.set_yticklabels([round(x,3) for x in yticks[2:8]])
#ax.set_ylim([-0.3,5])
#1/0
ax.xaxis.set_major_formatter(pter)

#Visual tweaks
ax.spines[['right', 'top']].set_visible(False)
for axis in ['bottom','left']:
    ax.spines[axis].set_linewidth(2)

ax.set_ylim(-0.3,1.08)
fig.tight_layout()
plt.xlabel(xlab_string + " (MB)")
plt.ylabel('LD')
plt.savefig('{0}.png'.format(savename), dpi=350,  bbox_inches='tight')
plt.savefig('{0}.svg'.format(savename), dpi=350,  bbox_inches='tight')
plt.show()
