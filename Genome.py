import multiprocessing.pool
from Bio import Entrez 
from Bio import SeqIO
import Templates
import multiprocessing
import time

BOX10 =  [['T', 0.8], ['A', 0.95], ['T', 0.45], ['A', 0.60], ['A', 0.50], ['T', 0.96]]
BOX35 =  [['T', 0.82], ['T', 0.84], ['G', 0.78], ['A', 0.65], ['C', 0.54], ['A', 0.45]]

REVERSE_DICT = {'A':'T', 'T':'A', 'G':'C', 'C':'G'}

AMINOACIDS = ["Ala", "Arg", "Asn", "Asp", "Cys", "Glu", "Gln", "Gly", "His", "Ile", "Leu", "Lys", "Met", "Phe", "Pro", "Ser", "Thr", "Trp", "Tyr", "Val"]


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
                start = time.time()
                
                with multiprocessing.Pool(4) as pool:
                    # Não roda 100% paralelo para evitar sobrecarga de requisições
                    results = [
                        pool.apply_async(self._GetrRNATemplate, ("5S",)),
                        pool.apply_async(self._GetrRNATemplate, ("16S",)),
                        pool.apply_async(self._GetrRNATemplate, ("23S",)),
                        pool.apply_async(self._GettRNAsTemplates)
                    ]

                    self.fiveSRNATemplate = results[0].get()
                    self.sixteenSRNATemplate = results[1].get()
                    self.twentyThreeSRNATemplate = results[2].get()
                    self.tRNAsTemplates = results[3].get()


                end = time.time()
                print(f"Tempo total: {end-start}")
                # self.sixteenSRNATemplate = self.__GetrRNATemplate("16S")
                # self.fiveSRNATemplate = self.__GetrRNATemplate("5S")
                # self.twentyThreeSRNATemplate = self.__GetrRNATemplate("23S")
                # self.tRNAsTemplates = self.__GettRNAsTemplates()

            self.tRNAsGenes = dict()    
            self.cutOff = rnaGenesCutoff
            self.rnaMinScore = rnaGenesMinScore

    def _GetrRNATemplate(self, RNAtype : str):
        Entrez.email = "dezinho_dh@hotmail.com"
        Entrez.api_key = "d13f56e78009d649f4f0dc96731dfaa8a308"
        orgnName = self.searchString if self.searchString.replace(" ", "").isalnum() else "Escherichia coli"
        upperLimit = 0; lowerLimit = 0
        match RNAtype:
            case "16S":              
                upperLimit = 1570; lowerLimit = 1530
            case "5S":
                upperLimit = 128; lowerLimit = 112
            case "23S":
                upperLimit = 2936; lowerLimit = 2876

        firstQueryString=f'''({RNAtype} rRNA[Title] AND {orgnName}[Orgn]) AND ("{lowerLimit-30}"[SLEN] : "{upperLimit+30}"[SLEN])'''
        secondQueryString = f"""(("{orgnName}"[Organism] OR {orgnName}[All Fields]) AND gene[All Fields] AND {RNAtype}[All Fields] AND rRNA[All Fields]) AND ("{lowerLimit-30}"[SLEN] : "{upperLimit+30}"[SLEN])"""

        template, firstRecords, secondRecords = self.__EntrezSearchrRNATemplate(firstQueryString, secondQueryString, upperLimit, lowerLimit)

        return template if template != "sem resultados" else self.__SearchWithHigherLimits(firstRecords, secondRecords, upperLimit+30, lowerLimit-30)
           
    def __EntrezSearchrRNATemplate(self, firstQueryString : str, secondQueryString : str, upperLimit : int, lowerLimit : int) -> str:
        try:
            handle = Entrez.esearch(db="nucleotide", term=firstQueryString, retmax=100)
            record = Entrez.read(handle)
            handle.close()
                
            handle = Entrez.efetch(db="nucleotide", id=record["IdList"], rettype="fasta")
            firstRecords = handle.read().split('>')
            handle.close()

            for result in firstRecords:
                result = result[result.find('\n'):].replace("\n", "")
                teste = len(result)
                if  upperLimit > len(result) > lowerLimit:
                    return result, None, None

            handle = Entrez.esearch(db="nucleotide", term=secondQueryString, retmax=100)
            record = Entrez.read(handle)
            handle.close()

            handle = Entrez.efetch(db="nucleotide", id=record["IdList"], rettype="fasta")
            secondRecords = handle.read().split('>')
            handle.close()

            for result in secondRecords:
                result = result[result.find('\n'):].replace("\n", "")
                teste = len(result)
                if  upperLimit > len(result) > lowerLimit:
                    return result, None, None
                
            return "sem resultados", firstRecords, secondRecords
        
        except Exception as e:
            print("Erro ao buscar pelos rRNAs do organismo selecionado, tente com outro ou inicialize sem searchString para utilizar os RNAs de Escherichia coli.\nMensagem de erro:", e)
        
    def __SearchWithHigherLimits(self, firstRecords : list, secondRecords : list, upperLimit : str, lowerLimit : str) -> str:
        for result in firstRecords:
            result = result[result.find('\n'):].replace('\n', "")
            teste = len(result)
            if upperLimit > len(result) > lowerLimit:
                return result
            
        for result in secondRecords:
            result = result[result.find('\n'):].replace('\n', "")
            teste = len(result)
            if upperLimit > len(result) > lowerLimit:
                return result
        return "sem resultados"     
        
        
    def _GettRNAsTemplates(self) -> dict:
        Entrez.email = "dezinho_dh@hotmail.com"
        Entrez.api_key = "d13f56e78009d649f4f0dc96731dfaa8a308"

        queryString = self.searchString if self.searchString.replace(" ", "").isalnum() else "Escherichia coli"
        dicttRNAs = {}

        try:
            baseGenome = self.__EntrezGetBaseGenome(queryString)

            for aminoacid in AMINOACIDS:
                result = self.__ExtracttRNASequence(baseGenome, aminoacid)

                teste = len(result)
                if (100 > len(result) > 65 ):
                    dicttRNAs[aminoacid] = result  
                else:
                    dicttRNAs[aminoacid] = "sem resultados"                   
                 
        except Exception as e:
            print("Erro ao buscar pelos tRNAs do organismo selecionado, tente com outro ou inicialize sem searchString para utilizar os RNAs de Escherichia coli.\nMensagem de erro:", e)
            
        return dicttRNAs

    def __EntrezGetBaseGenome(self, queryString) -> str:
        handle = Entrez.esearch(db="nucleotide", term=f"(tRNA[Feature key]) AND {queryString}[Title] AND (complete genome[Title] OR complete sequence[Title] NOT plasmid[Title] NOT fragment[Title] NOT partial[Title])", retmax=100)
        records = Entrez.read(handle)
        handle.close()

        details = self.__FetchDetails(records["IdList"])
        idsDates = [(result["Id"], result["CreateDate"]) for result in details]
        idsDates.sort(key= lambda idDate : idDate[1], reverse=True)

        handle = Entrez.efetch(db="nucleotide", id=idsDates[0][0], rettype="gb", retmode="text")
        record = SeqIO.read(handle, "genbank")
        handle.close()

        return record

    def __FetchDetails(self, ids):
        handle = Entrez.esummary(db="nucleotide", id=",".join(ids), retmode="xml")
        records = Entrez.read(handle)
        handle.close()
        return records

    def __ExtracttRNASequence(self, baseGenome, aminoacid):
        for feature in baseGenome.features:
            if feature.type == "tRNA" and feature.qualifiers["product"][0] == f"tRNA-{aminoacid}":
                start = feature.location.start
                end = feature.location.end
                return  str(baseGenome.seq[start:end])

        return "sem resultados"

    def SearchRNAGenes(self):
        start = time.time()
        self.__SearchrRNAGenes()
        self.__SearchtRNAGenes()
        end = time.time()
        print(f"Tempo de busca por todos os genes de rRNA: {end-start}")    
        return 
    
    def __SearchrRNAGenes(self):
        start = time.time()

        sixteenGenes = self._SinglerRNASearch(self.sixteenSRNATemplate)
        fiveGenes = self._SinglerRNASearch(self.fiveSRNATemplate)
        twentyThreeGenes = self._SinglerRNASearch(self.twentyThreeSRNATemplate)

        # with multiprocessing.Pool(3) as pool:
        #     results = [
        #     pool.apply_async(self._SinglerRNASearch, args= (self.sixteenSRNATemplate,)),
        #     pool.apply_async(self._SinglerRNASearch, args= (self.fiveSRNATemplate,)),
        #     pool.apply_async(self._SinglerRNASearch, args= (self.twentyThreeSRNATemplate,))
        #     ]

        #     sixteenGenes = results[0].get()
        #     fiveGenes = results[1].get()
        #     twentyThreeGenes = results[2].get()


        end = time.time()
        print(f"Tempo busca rRNA genes: {end-start}")

        sixteenGenes.sort(key=lambda gene : gene[2], reverse=True)
        fiveGenes.sort(key=lambda gene : gene[2], reverse=True)
        twentyThreeGenes.sort(key=lambda gene : gene[2], reverse=True)

        self.sixteenGenes = sixteenGenes
        self.fiveGenes = fiveGenes
        self.twentyThreeGenes = twentyThreeGenes

        return [sixteenGenes, fiveGenes, twentyThreeGenes]
    
    def _SinglerRNASearch(self, template : str) -> list:

        if not (template == "sem resultados"):       
            genes = self.__rRNAGenesIterations(template)
            if len(genes) == 0:
                genes = self.__rRNAGenesIterations(template.replace('N','A'))
        return genes

    
    def __rRNAGenesIterations(self, targetTemplate : str) -> list:
        number=0
        genes = list()

        while True:            
                gene = self.__GetRNAGene(number, targetTemplate)
                if isinstance(gene, str):
                    break
                number = number + self.cutOff
                if gene[2] > self.rnaMinScore:
                    genes.append(gene)
        return genes
        
    
    def __GetRNAGene(self, number : int, targetTemplate : str) -> list:
        sequence = self.sequence
        start = sequence.find(targetTemplate[:3], number)
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
        hasTemplatesList = list()

        teste = [(self.tRNAsTemplates[aminoacid], aminoacid) for aminoacid in AMINOACIDS]

        #Mudar para atribuir a resultado intermediário
        with multiprocessing.Pool(20) as pool:
            results = pool.starmap(self._SingletRNA, [(self.tRNAsTemplates[aminoacid], aminoacid) for aminoacid in AMINOACIDS], chunksize=20)
        
        i=0
        for aminoacid in AMINOACIDS:
            self.tRNAsGenes[aminoacid] = results[i]
            i += 1                         

        return
    
    def _SingletRNA(self, template: str, aminoacid : str) -> list:
        
        if template == "sem resultados":
            return "sem resultados"
        n = 0

        sequenceLength = len(self.sequence)
        currentTemplateLength = len(template)

        while True:            
                gene = self.__GetRNAGene(n, self.tRNAsTemplates[aminoacid])
                if isinstance(gene, str) or (n + self.cutOff) > (sequenceLength - currentTemplateLength):
                    break
                n = n + self.cutOff
                if gene[2] > self.rnaMinScore:
                    self.tRNAsGenes[aminoacid].append(gene)
        self.tRNAsGenes[aminoacid].sort(key= lambda gene : gene[2], reverse=True)  


    def __GetSequence(self, arch):
        arch = open(arch)
        arch.readline()
        sequence = ""
        for line in arch.readlines():
            sequence += line.strip()
        return sequence
    
    def __FindGenes(self, sequence : str, listGenes : list):
        n=0
        while True:       
            gene = self.__GetGene(sequence, n)
            if isinstance(gene, str): #Se for string, terminou a busca nessa sequência
                break   
            if gene[1]-gene[0] >= self.minimumSize: #Se não tiver o tamanho mínimo, pula esse gene
                geneBoxes = self.__GetBestBoxes(sequence, gene[0]) 
                product = geneBoxes[0]*geneBoxes[1]
                if(product >= self.minimumProbability):   
                    listGenes.append([gene[0], gene[1], geneBoxes[0], geneBoxes[1], product, geneBoxes[2]])   
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

        forwardStrandGenes = self.__FindGenes(sequence, forwardStrandGenes)
        reverseStrandGenes = self.__FindGenes(reverseSequence, reverseStrandGenes)
        
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
            print(f"{gene[0]}\t{gene[1]}\t{gene[4]}")
            print(f"Sequência do códon de início até o TATA box -35:")
            print(f"{gene[5]}")

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
# print(start-end)l