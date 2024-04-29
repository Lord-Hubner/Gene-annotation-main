import math
from Bio import Entrez 

BOX10 =  [['T', 0.8], ['A', 0.95], ['T', 0.45], ['A', 0.60], ['A', 0.50], ['T', 0.96]]
BOX35 =  [['T', 0.82], ['T', 0.84], ['G', 0.78], ['A', 0.65], ['C', 0.54], ['A', 0.45]]

REVERSE_DICT = {'A':'T', 'T':'A', 'G':'C', 'C':'G'}




class Genome:
    def __init__(self, arch : str, minProb = 0, minSize = 3 , includeRNAGenes = False, searchString = "") -> None:
        self.sequence = self.GetSequence(arch)
        self.minimumProbability = minProb
        self.minimumSize = minSize
        self.includeRNAGenes = includeRNAGenes
        if includeRNAGenes:
            self.searchString = searchString
            self.sixteenSRNA = self.SearchSequence("16S")
            self.fiveSRNA = self.SearchSequence("5S")
            self.twentythreeSRNA = self.SearchSequence("23S")

    def SearchSequence(self, RNAtype : str):
        queryString = self.searchString if self.searchString != "" else "E. coli"
        match RNAtype:
            case "16S":
                queryString=queryString+" 16S rRNA"
            case "5S":
                queryString=queryString+" 5S rRNA"
            case "23S":
                queryString=queryString+" 23S rRNA"
        self.queryString = queryString

        Entrez.email = "dezinho_dh@hotmail.com"
        try:
            handle = Entrez.esearch(db="nucleotide", term=queryString)
            record = Entrez.read(handle)
            handle.close()
        
            handle = Entrez.efetch(db="nucleotide", id=record["IdList"][0], rettype="fasta")
            record = handle.read()
            handle.close()


            return record
        except Exception as e:
            print("Erro ao buscar pelos rRNAs da sequência, mensagem de erro:", e)
            return None
        

    # GettRNAS(self, sequence : str):





    


    def GetSequence(self, arch):
        arch = open(arch)
        arch.readline()
        sequence = ""
        for line in arch.readlines():
            sequence += line.strip()
        return sequence


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
    
    # def GetRNAGenes(self, sequence: str):
    #     16sRNA
    
    # def Get16SRNAGene(self, ):
        

        

    def GetPontuation(self, boxSequence : str, boxProbs : list):
        points = 1
        i = 0
        for char in boxSequence:
            points *= boxProbs[i][1] if char == boxProbs[i][0] else (1-boxProbs[i][1])/3
            i = i + 1

        return points**(1/6)

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
            self.GetRNAGenes(sequence)
        return genes

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

    


# sequence = GetSequence("phageX174.fasta")
# result = GetStartEnd(sequence, 0)
# print(result)

# tata = 'TATAAT'
# sequence = GetSequence("phageX174.fasta")
# variable = GetBestBoxes(sequence, GetGene(sequence, 0)[0])
# print(variable[2])

genome = Genome("mycoplasmagenitalium.fasta", includeRNAGenes=True)
print(genome.AnnotateGenome())
print("bah men")