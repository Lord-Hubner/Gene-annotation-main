import math
from Bio import Entrez 
from numba import njit
from numba.experimental import jitclass
from numba import int32, string, boolean
import Templates
import time

BOX10 =  [['T', 0.8], ['A', 0.95], ['T', 0.45], ['A', 0.60], ['A', 0.50], ['T', 0.96]]
BOX35 =  [['T', 0.82], ['T', 0.84], ['G', 0.78], ['A', 0.65], ['C', 0.54], ['A', 0.45]]

REVERSE_DICT = {'A':'T', 'T':'A', 'G':'C', 'C':'G'}

AMINOACIDS = ["Ala", "Arg", "Asn", "Asp", "Cys", "Glu", "Gln", "Gly", "Hys", "Ile", "Leu", "Lys", "Met", "Phe", "Pro", "Ser", "Thr", "Trp", "Tyr", "Val"]


spec = [
    ('sequence', string),
    ('minimumProbability', int32),
    ('minimumSize', int32),
    ('includeRNAGenes', boolean)

]

#@jitclass
class Genome:
    #@njit
    def __init__(self, arch : str, minProb = 0, minSize = 3 , includeRNAGenes = False, searchString = "") -> None:
        '''
        Inicializa um genoma tendo como o alvo apontado por arch e busca por genes com o produto das probabildades dos promotores de pelo menos minProb
        e sequência de pelo menos minSize. Se includeRNAGenes for igual a true, busca também por genes de rRNA e tRNA no genoma tendo como base as sequências 
        do organismo específicado em searchString. Caso este último não seja específicado, a busca será em Escherichia coli.
        '''
        self.sequence = self.GetSequence(arch)
        self.minimumProbability = minProb
        self.minimumSize = minSize
        self.includeRNAGenes = includeRNAGenes
        if includeRNAGenes:
            self.searchString = searchString
            self.sixteenSRNA = self.GetrRNATemplate("16S")
            self.fiveSRNA = self.GetrRNATemplate("5S")
            self.twentythreeSRNA = self.GetrRNATemplate("23S")
            self.tRNAs = self.GettRNAsTemplates()

    #@njit
    def GetrRNATemplate(self, RNAtype : str):
        queryString = self.searchString if self.searchString != "" else "Escherichia coli"
        match RNAtype:
            case "16S":
                queryString=f"16S rRNA[Title] AND {queryString}[Orgn]"
            case "5S":
                queryString=f"5S rRNA[Title] AND {queryString}[Orgn]"
            case "23S":
                queryString=f"23S rRNA[Title] AND {queryString}[Orgn]"
        self.queryString = queryString

        Entrez.email = "dezinho_dh@hotmail.com"
        try:
            handle = Entrez.esearch(db="nucleotide", term=queryString, retmax=1)
            record = Entrez.read(handle)
            handle.close()
        
            handle = Entrez.efetch(db="nucleotide", id=record["IdList"][0], rettype="fasta")
            record = handle.read()
            handle.close()

            record = record[record.find('\n'):].replace("\n", "")
            return record
        except Exception as e:
            print("Erro ao buscar pelos rRNAs do organismo selecionado, tente com outro ou inicialize sem searchString para buscar em Escherichia coli.\nMensagem de erro:", e)
            return None
        
    #@njit
    def GettRNAsTemplates(self) -> dict:
        Entrez.email = "dezinho_dh@hotmail.com"
        queryString = self.searchString if self.searchString != "" else "Escherichia coli"
        dicttRNAs = {}

        for aminoacid in AMINOACIDS:
            records = self.EntrezSearchtRNA(aminoacid, queryString)

            try:
                for result in records[1:]:
                    result = result[result.find('\n'):].replace("\n", "")
                    if  (100 > len(result) > 65 ):                   
                        dicttRNAs[aminoacid] = result
                        break
            except Exception as e:
                print("Erro ao buscar pelos tRNAs do organismo selecionado, tente com outro ou inicialize sem searchString para buscar em Escherichia coli.\nMensagem de erro:", e)
                return None
            
        return dicttRNAs

    #@njit
    def EntrezSearchtRNA(self, aminoacid, queryString) -> str:
        try:
            handle = Entrez.esearch(db="nucleotide", term="tRNA-"+aminoacid+f"[Title] AND {queryString}[Orgn]", retmax=20)
            records = Entrez.read(handle)
            handle.close()

            handle = Entrez.efetch(db="nucleotide", id=records["IdList"], rettype="fasta")
            records = handle.read().split('>')
            handle.close()

            return records
        except Exception as e:
            print("Erro ao buscar pelos tRNAs do organismo selecionado, tente com outro ou inicialize sem searchString para buscar em Escherichia coli.\nMensagem de erro:", e)

    #@njit
    def SearchRNAGenes(self):
        sixteenrRNA, fiverRNA, twentythreeRNA = self.SearchrRNAGenes()

        tRNAs = self.SearchtRNAGenes()

        return [sixteenrRNA, fiverRNA, twentythreeRNA, tRNAs]
    
    #@njit
    def SearchtRNAGenes(self):
        dicttRNAs = dict

        for aminoacid in AMINOACIDS:
            dicttRNAs[aminoacid] = self.GetFivetRNAGenes(aminoacid)

    #@njit
    def GetFivetRNAGenes(self, aminoacid : str):
        sequence = self.sequence
        thisAminoacidTemplate = self.tRNAs[aminoacid][:3]

        genes = list()

        score = 0
        i = 0
        n = 0
        while i < 5:
            start = sequence.find(thisAminoacidTemplate, n) 

            for char in range(start, start+len(thisAminoacidTemplate)):
                if char == thisAminoacidTemplate[i]:
                    score = score + 1
           
            score = score/len(thisAminoacidTemplate)
            genes.append([start, char, score])
            n = char
            i = i + 1

        return genes

    #@njit
    def GetSRNAGenes(self, sequence, currentTargetTemplate):
        genes = list()
        n=0
        while True:
            gene = self.GetrRNAGene(sequence, n, currentTargetTemplate)
            if isinstance(gene, str):
                break
            n = gene[1]
            if gene[2] >= 0.75:
                genes.append(gene)

        return genes

    #@njit
    def GetrRNAGene(self, sequence : str, number : int, targetTemplate : str):

        start = sequence.find(targetTemplate[:3], number)
        end = start+len(targetTemplate)
        score = 0
        i=0

        for n in range(start, end):
            if sequence[n] == targetTemplate[i]:
                score = score + 1
            i = i + 1
        score = score/len(targetTemplate)
        return [start, n, score]

    #@njit
    def GetSequence(self, arch):
        arch = open(arch)
        arch.readline()
        sequence = ""
        for line in arch.readlines():
            sequence += line.strip()
        return sequence

    #@njit
    def GetGene(self, sequence: str, number : int):

        startIndex = sequence.find("ATG", number)
        if(startIndex == -1):
            return f"Nenhum possível gene a partir do nucleotídeo na posição {number}"
        for n in range(startIndex, len(sequence), 3):
            if sequence[n:n+3] == "TAG" or sequence[n:n+3] == "TAA" or sequence[n:n+3] == "TGA":
                stopIndex = n
                return [startIndex, stopIndex]
        #Pode buscar de novo começando desde o início já que o genoma é circular
        for n in range(0, len(sequence), 3):
            if sequence[n:n+3] == "TAG" or sequence[n:n+3] == "TAA" or sequence[n:n+3] == "TGA":
                stopIndex = n
                return [startIndex, stopIndex]

        return f"Nenhum possível gene a partir do nucleotídeo na posição {number}"

    #@njit
    def GetPontuation(self, boxSequence : str, boxProbs : list):
        points = 1
        i = 0
        for char in boxSequence:
            points *= boxProbs[i][1] if char == boxProbs[i][0] else (1-boxProbs[i][1])/3
            i = i + 1

        return points**(1/6)

    #@njit
    def GetBestBoxes(self, sequence : str, startCodonIndex : int):
        potential10Boxes = list()
        potential35Boxes = list()
        for x in range(startCodonIndex-15, startCodonIndex-11, 1):
            potential10Boxes.append([x, sequence[x:x+6]]) 
        
        listPoints = list()
        for box in potential10Boxes:
            listPoints.append([box[0], box[1], self.GetPontuation(box[1], BOX10)])

        best10Box = max(listPoints, key= lambda x: x[2])

        for x in range(best10Box[0]-25, best10Box[0]-22, 1):
            potential35Boxes.append([x, sequence[x:x+6]])
        
        listPoints.clear()
        for box in potential35Boxes:
            listPoints.append([box[0], box[1], self.GetPontuation(box[1], BOX35)])

        best35Box = max(listPoints, key= lambda x: x[2])

        start = sequence[best35Box[0]:startCodonIndex+3] if best35Box[0] >= 0 else sequence[best35Box[0]:]+sequence[:startCodonIndex+3]
        returnString = start+"\n"+"".join(a[0] for a in BOX35)+(best10Box[0]-(best35Box[0]+len(best35Box[1])))*" "+"".join(a[0] for a in BOX10)+(startCodonIndex-(best10Box[0]+len(best10Box[1])))*" "+"ATG"

        print(returnString)
        return [best10Box[2], best35Box[2], returnString]

    #@njit
    def AnnotateGenome(self):
        
        '''
        Dado o nome de um arquivo fasta na pasta local, uma probabilidade e um inteiro, retorna os possíveis genes
        da sequência desse arquivo que tenham o produto das probabilidades dos boxes de pelo menos cutProb
        e tamanho mínimo de minSize
        '''
        sequence = self.sequence
        translationTable = str.maketrans(REVERSE_DICT)
        reverseSequence = sequence.translate(translationTable)
        forwardStrandGenes = list()
        reverseStrandGenes = list()
        n = 0

        forwardStrandGenes = self.FindGenes(sequence, forwardStrandGenes, n)
        reverseStrandGenes = self.FindGenes(reverseSequence, reverseStrandGenes, n)
        
        genes = forwardStrandGenes
        genes.extend(reverseStrandGenes)
        genes.sort(key= lambda gene: gene[2], reverse=True)

        if(self.includeRNAGenes):
            sixteenSRNA, fiveSRNA, twentythreeSRNA, tRNAs = self.SearchRNAGenes()
            
        return [genes, sixteenSRNA, fiveSRNA, twentythreeSRNA, tRNAs]

    #@njit
    def FindGenes(self, sequence : str, listGenes : list, n : int):
        while True:       
            gene = self.GetGene(sequence, n)
            if isinstance(gene, str): #Se for string, terminou a busca nessa sequência
                break   
            if gene[1]-gene[0] >= self.minimumSize: #Se não tiver o tamanho mínimo, pula esse gene
                geneBoxes = self.GetBestBoxes(sequence, gene[0]) 
                product = geneBoxes[0]*geneBoxes[1]
                if(product >= self.minimumProbability):   
                    listGenes.append([gene[0], gene[1], product])   
            n = gene[0]+1
        return listGenes

start = time.time()
genome = Genome("mycoplasmagenitalium.fasta", includeRNAGenes=True)
print(genome.AnnotateGenome())
print("bah men")
end = time.time()
print(start-end)