from genome import Genome

noFilters = Genome("phageX174.fasta")


minPercent10 = Genome("phageX174.fasta", 0.10)

minPercent30 = Genome("phageX174.fasta", 0.30)

minPercent5 = Genome("phageX174.fasta", 0.05)

minSize30 = Genome("phageX174.fasta", 0, 30)

minSize100 = Genome("phageX174.fasta", 0, 100)

minSize300 = Genome("phageX174.fasta", 0, 300)

minPercent10minSize30 = Genome("phageX174.fasta", 0.10, 30)

minPercent5minSize30 = Genome("phageX174.fasta", 0.05, 30)

minPercent5minSize100 = Genome("phageX174.fasta", 0.05, 100)