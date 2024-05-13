from Bio import Entrez 
import Templates
import time

BOX10 =  [['T', 0.8], ['A', 0.95], ['T', 0.45], ['A', 0.60], ['A', 0.50], ['T', 0.96]]
BOX35 =  [['T', 0.82], ['T', 0.84], ['G', 0.78], ['A', 0.65], ['C', 0.54], ['A', 0.45]]

REVERSE_DICT = {'A':'T', 'T':'A', 'G':'C', 'C':'G'}

AMINOACIDS = ["Ala", "Arg", "Asn", "Asp", "Cys", "Glu", "Gln", "Gly", "Hys", "Ile", "Leu", "Lys", "Met", "Phe", "Pro", "Ser", "Thr", "Trp", "Tyr", "Val"]


class Genome:
    Initialized = False
    def __init__(self, arch : str, minProb = 0, minSize = 3 , includeRNAGenes = False, searchString = "", rnaGenesCutoff = 30, rnaGenesMinScore = 0) -> None:
        '''
        Inicializa um genoma tendo como o alvo apontado por arch e busca por genes com o produto das probabilidades dos promotores de pelo menos minProb
        e sequência de tamanho de pelo menos minSize. Se includeRNAGenes for igual a true, busca também por genes de rRNA e tRNA no genoma tendo como base 
        as sequências do organismo específicado em searchString e permitindo identificar genes a rnaGenesCutoff nucleotídeos de distância. Caso searchString
        não seja específicado, serão utilizadas sequências locais de Escherichia coli. Caso seja, será utilizado o email especificado em email para fazer as 
        consultas no NCBI


        Parâmetros:
        arch (str): Caminho do arquivo fasta contendo o genoma alvo.
        minProb (float, opcional): Produto mínimo das probabilidades dos promotores. Padrão 0.
        minSize (int, opcional): Tamanho mínimo do gene. Padrão 3.
        includeRNAGenes (bool, opcional): Se True, busca por genes de RNA também. Padrão False.
        searchString (str, opcional): Organismo cujos genes de RNA serão usados como alvo para a busca na sequência fornecida. Por padrão utiliza sequências locais de Escherichia coli.
        rnaGenesCutoff (int, opcional): Intervalo máximo entre genes para buscar genes de RNA. Valores muito baixos podem tornar a busca bem mais demorada e devolver muitos resultados repetidos. Padrão 30.
        '''

        self.sequence = self.__GetSequence(arch)
        self.minimumProbability = minProb
        self.minimumSize = minSize
        self.includeRNAGenes = includeRNAGenes
        if includeRNAGenes:
            self.searchString = searchString
            if searchString == "":
                self.sixteenSRNATemplate = Templates.sixteenSRNA 
                self.fiveSRNATemplate = Templates.fiveSRNA
                self.twentyThreeSRNATemplate = Templates.twentyThreeSRNA
                self.tRNAsTemplates = Templates.tRNAs
            else:
                self.sixteenSRNATemplate = self.__GetrRNATemplate("16S")
                self.fiveSRNATemplate = self.__GetrRNATemplate("5S")
                self.twentyThreeSRNATemplate = self.__GetrRNATemplate("23S")
                self.tRNAsTemplates = self.__GettRNAsTemplates()
            self.tRNAsGenes = dict()    
            self.cutOff = rnaGenesCutoff
            self.rnaMinScore = rnaGenesMinScore

    def __GetrRNATemplate(self, RNAtype : str):
        queryString = self.searchString if self.searchString.isalnum() else "Escherichia coli"
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
            print("Erro ao buscar pelos rRNAs do organismo selecionado, tente com outro ou inicialize sem searchString para utilizar os RNAs de Escherichia coli.\nMensagem de erro:", e)
            return None
        
    def __GettRNAsTemplates(self) -> dict:
        Entrez.email = "dezinho_dh@hotmail.com"
        queryString = self.searchString if self.searchString.isalnum() else "Escherichia coli"
        dicttRNAs = {}

        for aminoacid in AMINOACIDS:
            records = self.__EntrezSearchtRNA(aminoacid, queryString)

            try:
                for result in records[1:]:
                    result = result[result.find('\n'):].replace("\n", "")
                    if  (100 > len(result) > 65 ):                   
                        dicttRNAs[aminoacid] = result
                        break
            except Exception as e:
                print("Erro ao buscar pelos tRNAs do organismo selecionado, tente com outro ou inicialize sem searchString para utilizar os RNAs de Escherichia coli.\nMensagem de erro:", e)
                return None
            
        return dicttRNAs

    def __EntrezSearchtRNA(self, aminoacid, queryString) -> str:
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

    def SearchRNAGenes(self):
        self.__SearchtRNAGenes()
        self.__GetrRNAGenes()
        return 
    
    def __GetrRNAGenes(self):
        sixteenGenes = list()
        fiveGenes = list()
        twentyThreeGenes = list()

        n=0

        initiate = time.time()
        while True:
            gene = self.__GetRNAGene(n, self.sixteenSRNATemplate)
            if isinstance(gene, str) :
                break
            n =  n + self.cutOff 
            if gene[2] > self.rnaMinScore:
                sixteenGenes.append(gene)
        finalize = time.time()
        print(f"Time 16: {finalize - initiate}\n")
        
        n=0
        initiate = time.time()
        while True:            
            gene = self.__GetRNAGene(n, self.fiveSRNATemplate)
            if isinstance(gene, str) :
                break
            n = n + self.cutOff
            if gene[2] > self.rnaMinScore:
                fiveGenes.append(gene)
        finalize = time.time()
        print(f"Time 5: {finalize - initiate}\n")

        n=0
        initiate = time.time()
        while True:
            gene = self.__GetRNAGene(n, self.twentyThreeSRNATemplate)
            if isinstance(gene, str) :
                break
            n = n + self.cutOff
            if gene[2] > self.rnaMinScore:
                twentyThreeGenes.append(gene)
        finalize = time.time()
        print(f"Time 23: {finalize - initiate}\n")
        
        initiate = time.time()

        sixteenGenes.sort(key=lambda gene : gene[2], reverse=True)
        fiveGenes.sort(key=lambda gene : gene[2], reverse=True)
        twentyThreeGenes.sort(key=lambda gene : gene[2], reverse=True)

        self.sixteenGenes = sixteenGenes
        self.fiveGenes = fiveGenes
        self.twentyThreeGenes = twentyThreeGenes

        finalize = time.time()
        print(f"Time srt: {finalize - initiate}\n")

        return [sixteenGenes, fiveGenes, twentyThreeGenes]
    
    def __GetRNAGene(self, number : int, targetTemplate : str):
        sequence = self.sequence
        start = sequence.find(targetTemplate[:3].replace('N', 'A'), number)
        if(start == -1):
            return f"Nenhum possível gene de rRNA a partir do nucleotídeo na posição {number}"
        end = start+len(targetTemplate)
        score = 0
        i=0

        sequence = sequence+sequence[:len(targetTemplate)+1]

        for n in range(start, end):
            if sequence[n] == targetTemplate[i]:
                score = score + 1
            i = i + 1
        score = score/len(targetTemplate)
        return [start, n, score]
    
    def __SearchtRNAGenes(self):
        sequenceLength = len(self.sequence)

        initiate = time.time()
        for aminoacid in AMINOACIDS:
            currentTemplateLength = len(self.tRNAsTemplates[aminoacid])
            self.tRNAsGenes[aminoacid] = list()
            n = 0
            while True:            
                 gene = self.__GetRNAGene(n, self.tRNAsTemplates[aminoacid])
                 if isinstance(gene, str) or (n + self.cutOff) > (sequenceLength - currentTemplateLength):
                    break
                 n = n + self.cutOff
                 if gene[2] > self.rnaMinScore:
                    self.tRNAsGenes[aminoacid].append(gene)
            self.tRNAsGenes[aminoacid].sort(key= lambda gene : gene[2], reverse=True)
                
        finalize = time.time()
        print(f"Time tRNAs: {finalize - initiate}\n")

        return

    def __GetSequence(self, arch):
        arch = open(arch)
        arch.readline()
        sequence = ""
        for line in arch.readlines():
            sequence += line.strip()
        return sequence
    
    def __FindGenes(self, sequence : str, listGenes : list, n : int):
        while True:       
            gene = self.__GetGene(sequence, n)
            if isinstance(gene, str): #Se for string, terminou a busca nessa sequência
                break   
            if gene[1]-gene[0] >= self.minimumSize: #Se não tiver o tamanho mínimo, pula esse gene
                geneBoxes = self.__GetBestBoxes(sequence, gene[0]) 
                product = geneBoxes[0]*geneBoxes[1]
                if(product >= self.minimumProbability):   
                    listGenes.append([gene[0], gene[1], geneBoxes[0], geneBoxes[1], product])   
            n = gene[0]+1
        return listGenes
    
    def __GetGene(self, sequence: str, number : int):

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
    
    def __GetBestBoxes(self, sequence : str, startCodonIndex : int):
        potential10Boxes = list()
        potential35Boxes = list()
        for x in range(startCodonIndex-15, startCodonIndex-11, 1):
            potential10Boxes.append([x, sequence[x:x+6]]) 
        
        listPoints = list()
        for box in potential10Boxes:
            listPoints.append([box[0], box[1], self.__GetPontuation(box[1], BOX10)])

        best10Box = max(listPoints, key= lambda x: x[2])

        for x in range(best10Box[0]-25, best10Box[0]-22, 1):
            potential35Boxes.append([x, sequence[x:x+6]])
        
        listPoints.clear()
        for box in potential35Boxes:
            listPoints.append([box[0], box[1], self.__GetPontuation(box[1], BOX35)])

        best35Box = max(listPoints, key= lambda x: x[2])

        start = sequence[best35Box[0]:startCodonIndex+3] if best35Box[0] >= 0 else sequence[best35Box[0]:]+sequence[:startCodonIndex+3]
        returnString = start+"\n"+"".join(a[0] for a in BOX35)+(best10Box[0]-(best35Box[0]+len(best35Box[1])))*" "+"".join(a[0] for a in BOX10)+(startCodonIndex-(best10Box[0]+len(best10Box[1])))*" "+"ATG"

        # print(returnString)
        return [best10Box[2], best35Box[2], returnString]

    def __GetPontuation(self, boxSequence : str, boxProbs : list):
        points = 1
        i = 0
        for char in boxSequence:
            points *= boxProbs[i][1] if char == boxProbs[i][0] else (1-boxProbs[i][1])/3
            i = i + 1

        return points**(1/6)

    def AnnotateGenome(self):
        
        '''
        Roda a anotação do Genoma de acordo com os parâmetros específicados no construtor
        '''
        sequence = self.sequence
        translationTable = str.maketrans(REVERSE_DICT)
        reverseSequence = sequence.translate(translationTable)
        forwardStrandGenes = list()
        reverseStrandGenes = list()
        n = 0

        forwardStrandGenes = self.__FindGenes(sequence, forwardStrandGenes, n)
        reverseStrandGenes = self.__FindGenes(reverseSequence, reverseStrandGenes, n)
        
        genes = forwardStrandGenes
        genes.extend(reverseStrandGenes)
        genes.sort(key= lambda gene: gene[4], reverse=True)

        if(self.includeRNAGenes):
            self.SearchRNAGenes()

        self.DNAGenes = genes
        self.Initialized = True
       
        if self.includeRNAGenes:
            return [genes, self.sixteenGenes, self.fiveGenes, self.twentyThreeGenes, self.tRNAsGenes]
        
        return genes

    def PrintResults(self):
        '''
        Retorna uma representação visual dos resultados da última anotação de genoma realizada.
        '''
        if not self.Initialized:
            raise Exception("Rode AnnotateGenomes() antes!")

        print("Possíveis genes de DNA:")
        print("Início\tFim\tProduto da pontuação")
        for gene in self.DNAGenes:
            print(f"{gene[0]}\t{gene[1]}\t{gene[2]}")

        if self.includeRNAGenes:
            print()
            print("Possíveis genes de 16S rRNA:")
            print("Início\tFim\tScore")
            for gene in self.sixteenGenes:
                print(f"{gene[0]}\t{gene[1]}\t{gene[2]}")
            
            print()
            print("Possíveis genes de 5S rRNA:")
            print("Início\tFim\tScore")
            for gene in self.fiveGenes:
                print(f"{gene[0]}\t{gene[1]}\t{gene[2]}")

            print()
            print("Possíveis genes de 23S rRNA:")
            print("Início\tFim\tScore")
            for gene in self.twentyThreeGenes:
                print(f"{gene[0]}\t{gene[1]}\t{gene[2]}")
        
            for aminoacid in AMINOACIDS:
                print()
                print(f"Possíveis genes de tRNA-{aminoacid}:")
                print("Início\tFim\tScore")
                for gene in self.tRNAsGenes[aminoacid]:
                    print(f"{gene[0]}\t{gene[1]}\t{gene[2]}")


# start = time.time()
# genome = Genome("mycoplasmagenitalium.fasta", includeRNAGenes=True)
# genome.AnnotateGenome()
# genome.PrintResults()
# end = time.time()
# print(start-end)