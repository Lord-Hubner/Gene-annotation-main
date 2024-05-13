from Genome import Genome
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
import math

noFilters = Genome("mycoplasmagenitalium.fasta")
noFilters.AnnotateGenome()

minPercent5 = Genome("mycoplasmagenitalium.fasta",  0.05)
minPercent5.AnnotateGenome()

minPercent10 = Genome("mycoplasmagenitalium.fasta",  0.10)
minPercent10.AnnotateGenome()

minPercent15 = Genome("mycoplasmagenitalium.fasta",  0.15)
minPercent15.AnnotateGenome()

minPercent20 = Genome("mycoplasmagenitalium.fasta",  0.20)
minPercent20.AnnotateGenome()

minPercent25 = Genome("mycoplasmagenitalium.fasta",  0.25)
minPercent25.AnnotateGenome()

minPercent30 = Genome("mycoplasmagenitalium.fasta", 0.30)
minPercent30.AnnotateGenome()


nGenesNoFilters = len(noFilters.DNAGenes)
nGenesPercent5 =  len(minPercent5.DNAGenes)
nGenesPercent10 = len(minPercent10.DNAGenes)
nGenesPercent15 =  len(minPercent15.DNAGenes)
nGenesPercent20 =  len(minPercent20.DNAGenes)
nGenesPercent25 =  len(minPercent25.DNAGenes)
nGenesPercent30 = len(minPercent30.DNAGenes)

#Número de genes encontrados por tamanho mínimo de CDS

minSize30 = Genome("mycoplasmagenitalium.fasta", 0, 30)
minSize30.AnnotateGenome()

minSize100 = Genome("mycoplasmagenitalium.fasta", 0, 100)
minSize100.AnnotateGenome()

minSize300 = Genome("mycoplasmagenitalium.fasta", 0, 300)
minSize300.AnnotateGenome()

minSize600 = Genome("mycoplasmagenitalium.fasta", 0, 600)
minSize600.AnnotateGenome()

minSize900 = Genome("mycoplasmagenitalium.fasta", 0, 900)
minSize900.AnnotateGenome()

minSize1200 = Genome("mycoplasmagenitalium.fasta", 0, 1200)
minSize1200.AnnotateGenome()

minSize1500 = Genome("mycoplasmagenitalium.fasta", 0, 1500)
minSize1500.AnnotateGenome()

nGenesminSize30 = len(minSize30.DNAGenes)
nGenesminSize100 = len(minSize100.DNAGenes)
nGenesminSize300 = len(minSize300.DNAGenes)
nGenesminSize600 = len(minSize600.DNAGenes)
nGenesminSize900 = len(minSize900.DNAGenes)
nGenesminSize1200 = len(minSize1200.DNAGenes)
nGenesminSize1500 = len(minSize1500.DNAGenes)

y = [nGenesNoFilters, nGenesminSize30, nGenesminSize100, nGenesminSize300, nGenesminSize600, nGenesminSize900, nGenesminSize1200, nGenesminSize1500]
x = [0, 30, 100] + list(np.arange(300, 1800, 300))

fig, (ax1, ax2) = plt.subplots(2)
fig.suptitle("Diferença de número de genes obtidos por filtro de tamanho\nMycoplasma genitalium")
fig.set_size_inches([8, 8])
ax1.bar(x, y, width=20)
ax1.set_xlabel("Tamanho mínimo")
ax1.set_xticks(x)
ax1.set_xticklabels(x)
ax1.set_ylabel("Número de genes")

ax2.plot(x, y, '-o')
ax2.set_yscale("log")
ax2.set_xlabel("Tamanho mínimo")
ax2.set_ylabel("log(número de genes)")
ax2.set_xticks(x)
ax2.set_yticks(y)
ax2.set_xticklabels(x)
ax2.set_yticklabels(y)
plt.show()

###Probabilidades BOX10 e BOX35 dos genes obtidos com produto mínimo 0 até 30

y = [nGenesNoFilters, nGenesPercent5, nGenesPercent10, nGenesPercent15, nGenesPercent20, nGenesPercent25, nGenesPercent30]

fig, axes = plt.subplots(3,2)

genes = [a[2] for a in noFilters.DNAGenes]
axes[0][0].hist(genes, 50, color='darkblue')
# plt.xticks(x, ["Sem filtro"]+ [a for a in np.arange(5, 35, 5)])
# axes[0][0].set_ylabel("Número de genes")
axes[0][0].set_title("BOX10 - sem filtro")

genes = [a[3] for a in noFilters.DNAGenes]
axes[0][1].hist(genes, 50, color='darkblue')
# plt.xticks(x, ["Sem filtro"]+ [a for a in np.arange(5, 35, 5)])
axes[0][1].set_title("BOX35 - sem filtro")

genes = [a[2] for a in minPercent15.DNAGenes]
axes[1][0].hist(genes, 50, color='yellow')
# plt.xticks(x, ["Sem filtro"]+ [a for a in np.arange(5, 35, 5)])
# axes[1][0].set_ylabel("Número de genes")
axes[1][0].set_title("BOX10 - produto mínimo de 15%")

genes = [a[3] for a in minPercent15.DNAGenes]
axes[1][1].hist(genes, 50, color='yellow')
# plt.xticks(x, ["Sem filtro"]+ [a for a in np.arange(5, 35, 5)])
axes[1][1].set_title("BOX35 - produto mínimo de 15%")

genes = [a[2] for a in minSize600.DNAGenes]
axes[2][0].hist(genes, 50, color='olive')
# plt.xticks(x, ["Sem filtro"]+ [a for a in np.arange(5, 35, 5)])
# axes[2][0].set_xlabel("Probabilidade")
# axes[2][0].set_ylabel("Número de genes")
axes[2][0].set_title("BOX10 - tamanho mínimo de 600bp")

genes = [a[3] for a in minSize600.DNAGenes]
axes[2][1].hist(genes, 50, color='olive')
# plt.xticks(x, ["Sem filtro"]+ [a for a in np.arange(5, 35, 5)])
# axes[2][1].set_xlabel("Probabilidade")
axes[2][1].set_title("BOX35 - tamanho mínimo de 600bp")

fig.set_size_inches([10, 8])
fig.suptitle("Número de genes por probabilidade do Box")
fig.text(0.5, 0.02, 'Probabilidade', ha='center')
fig.text(0.02, 0.5, 'Número de genes', va='center', rotation='vertical')
plt.subplots_adjust(hspace=0.4)
plt.show()

###Tamanho e número das possíveis proteínas identificadas 

fig, axes = plt.subplots(1,3)

genes = [a[1]-a[0] for a in noFilters.DNAGenes]
axes[0].hist(genes, 50, color='forestgreen')
axes[0].set_title("Sem filtro")

genes = [a[1]-a[0] for a in minSize600.DNAGenes]
axes[1].hist(genes, 50, color='forestgreen')
axes[1].set_title("Tamanho mínimo 600")

genes = [a[1]-a[0] for a in minSize1500.DNAGenes]
axes[2].hist(genes, 50, color='forestgreen')
axes[2].set_title("Tamanho mínimo 1500")
axes[2].yaxis.set_major_locator(ticker.MultipleLocator(1))

fig.suptitle("Número de genes por tamanho")

fig.text(0.5, 0.02, 'Tamanho do gene', ha='center')
fig.text(0.02, 0.5, 'Número de genes', va='center', rotation='vertical')
plt.show()

#Tamanho mínimo de 300 e produto de 15%
minPercent15Size300 = Genome("mycoplasmagenitalium.fasta", 0.15, 300)
minPercent15Size300.AnnotateGenome()

genes = [a[4] for a in minPercent15Size300.DNAGenes]
plt.title("Genes que são bons candidatos")
plt.xlabel("Produto das probabilidades")
plt.ylabel("Número de genes")
axes = plt.gca()
axes.yaxis.set_major_locator(ticker.MultipleLocator(1))

plt.hist(genes, 50, color='red')
plt.show()


#RNAs
#Números de genes de rRNA com o cutoff padrão vs número de genes de rRNA com score > 0.25

rnasDefault = Genome("mycoplasmagenitalium.fasta", includeRNAGenes=True)
rnasDefault.SearchRNAGenes()

rnasCutoff = Genome("mycoplasmagenitalium.fasta", includeRNAGenes=True, rnaGenesMinScore=0.25)
rnasCutoff.SearchRNAGenes()

fig, axes = plt.subplots(3,2)

genes = [a[2] for a in rnasDefault.sixteenGenes]
axes[0][0].hist(genes, 50, color='yellowgreen')
axes[0][0].set_title("Sem filtro")
axes[0][0].set_ylabel("16S rRNA")

genes = [a[2] for a in rnasDefault.fiveGenes]
axes[1][0].hist(genes, 50, color='yellowgreen')
axes[1][0].set_ylabel("5S rRNA")

genes = [a[2] for a in rnasDefault.twentyThreeGenes]
axes[2][0].hist(genes, 50, color='yellowgreen')
axes[2][0].set_ylabel("23S rRNA")


genes = [a[2] for a in rnasCutoff.sixteenGenes]
axes[0][1].hist(genes, 50, color='darkgreen')
axes[0][1].set_title("Score > 25%")

genes = [a[2] for a in rnasCutoff.fiveGenes]
axes[1][1].hist(genes, 50, color='darkgreen')

genes = [a[2] for a in rnasCutoff.twentyThreeGenes]
axes[2][1].hist(genes, 50, color='darkgreen')

fig.text(0.5, 0.02, 'Score do gene', ha='center')
fig.text(0.02, 0.5, 'Número de genes', va='center', rotation='vertical')
fig.suptitle("Número de possíveis genes de cada rRNA")
plt.show()

#Números de genes de tRNA com o cutoff padrão vs número de genes de tRNA com score > 0.25
#Aqui foram usados os tRNA-Ala, tRNA-Glu e tRNA-Thr como demonstrativos dos 20 resultados para tRNAs

fig, axes = plt.subplots(3,2)

genes = [a[2] for a in rnasDefault.tRNAsGenes["Ala"]]
axes[0][0].hist(genes, 50, color='yellowgreen')
axes[0][0].set_title("Sem filtro")
axes[0][0].set_ylabel("tRNA-Ala")

genes = [a[2] for a in rnasDefault.tRNAsGenes["Glu"]]
axes[1][0].hist(genes, 50, color='yellowgreen')
axes[1][0].set_ylabel("tRNA-Glu")

genes = [a[2] for a in rnasDefault.tRNAsGenes["Thr"]]
axes[2][0].hist(genes, 50, color='yellowgreen')
axes[2][0].set_ylabel("tRNA-Thr")


genes = [a[2] for a in rnasCutoff.tRNAsGenes["Ala"]]
axes[0][1].hist(genes, 50, color='darkgreen')
axes[0][1].set_title("Score > 25%")

genes = [a[2] for a in rnasCutoff.tRNAsGenes["Glu"]]
axes[1][1].hist(genes, 50, color='darkgreen')

genes = [a[2] for a in rnasCutoff.tRNAsGenes["Thr"]]
axes[2][1].hist(genes, 50, color='darkgreen')

fig.text(0.5, 0.02, 'Score do gene', ha='center')
fig.text(0.02, 0.5, 'Número de genes', va='center', rotation='vertical')
fig.suptitle("Número de possíveis genes de três tRNAs")
plt.show()


# minPercent10minSize30 = Genome("phageX174.fasta", 0.10, 30)

# minPercent5minSize30 = Genome("phageX174.fasta", 0.05, 30)

# minPercent5minSize100 = Genome("phageX174.fasta", 0.05, 100)