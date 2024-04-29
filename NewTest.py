from Bio import Entrez

# Provide your email address to NCBI
Entrez.email = "dezinho_dh@hotmail.com"

def search_16s_sequence(query):
    try:
        handle = Entrez.esearch(db="nucleotide", term=query)
        record = Entrez.read(handle)
        handle.close()
        return record["IdList"]
    except Exception as e:
        print("Error searching for sequence:", e)
        return None

def fetch_sequence_by_id(id):
    try:
        handle = Entrez.efetch(db="nucleotide", id=id, rettype="fasta")
        record = handle.read()
        handle.close()
        return record
    except Exception as e:
        print("Error fetching sequence:", e)
        return None

# Example: Search query for a specific organism name
organism_name = "rrna"

# Search for 16S rRNA sequences
search_result_ids = search_16s_sequence(organism_name)

for i in search_result_ids:
    # Fetch the first sequence from the search results
    sequence = fetch_sequence_by_id(i)
    
    if sequence:
        print("16S rRNA Sequence:")
        print(sequence)
    else:
        print("Sequence retrieval failed.")
else:
    print("No sequences found for the organism:", organism_name)