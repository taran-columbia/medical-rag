import os
from Bio import Entrez
import json
from dotenv import load_dotenv

load_dotenv()
# NCBI requires an email address to use their API
Entrez.email = os.getenv("EMAIL_USER", "default@example.com")

def fetch_pubmed_abstracts(query, max_results=20):
    print(f"Searching PubMed for: {query}...")
    
    # 1. Search for the query to get a list of IDs (PMIDs)
    search_handle = Entrez.esearch(db="pubmed", term=query, retmax=max_results)
    search_results = Entrez.read(search_handle)
    search_handle.close()
    
    id_list = search_results["IdList"]
    if not id_list:
        print("No results found.")
        return []

    # 2. Fetch the actual abstracts using the IDs
    fetch_handle = Entrez.efetch(db="pubmed", id=id_list, rettype="medline", retmode="text")
    from Bio import Medline
    records = Medline.parse(fetch_handle)
    
    abstracts = []
    for record in records:
        # We only want entries that actually have an abstract (AB) and a title (TI)
        if "AB" in record and "TI" in record:
            full_text = f"Title: {record['TI']}\nAbstract: {record['AB']}"
            abstracts.append(full_text)
            
    fetch_handle.close()
    print(f"✅ Successfully fetched {len(abstracts)} clinical abstracts.")
    return abstracts

if __name__ == "__main__":
    # Example: Fetch 10 real papers about Diabetes and Retinopathy
    data = fetch_pubmed_abstracts("Diabetes and Retinopathy", max_results=10)
    
    # Save to a local file so you don't have to call the API every time
    with open("raw_medical_data.json", "w") as f:
        json.dump(data, f)