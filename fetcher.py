import os
from Bio import Entrez
from Bio import Medline
from dotenv import load_dotenv
import re
load_dotenv()
Entrez.email = os.getenv("EMAIL_USER", "default@example.com") # Required by NCBI

def convert_phrases_to_plus_string(text):
    """
    Splits a string by commas and whitespace, then joins 
    the resulting words with a '+' sign.
    """
    # Split by one or more occurrences of either a comma or whitespace
    # This automatically handles "phrase 1, phrase 2" and "phrase1,phrase2"
    words = re.split(r'[,\s]+', text)
    
    # Filter out empty strings that might occur from trailing/leading separators
    clean_words = [word for word in words if word]
    
    # Join the final list with a '+'
    return '+'.join(clean_words)

def fetch_medical_papers(user_question, max_results=3):
    # Convert sentence to keyword string
    pubmed_query = convert_phrases_to_plus_string(user_question)
    print(f"  🔍 PubMed Query: {pubmed_query}")

    try:
        # 1. Search for IDs
        search_handle = Entrez.esearch(db="pubmed", term=pubmed_query, retmax=max_results)
        search_results = Entrez.read(search_handle)
        search_handle.close()

        id_list = search_results.get("IdList", [])
        print(f"  ✅ Found {len(id_list)} PMIDs")

        if not id_list:
            return []

        # 2. Fetch the actual content
        fetch_handle = Entrez.efetch(db="pubmed", id=id_list, rettype="medline", retmode="text")
        records = Medline.parse(fetch_handle)
        
        abstracts = []
        for r in records:
            title = r.get("TI", "No Title")
            # Fallback: if no abstract, use the title as context
            abstract = r.get("AB", f"Note: Full abstract not found. Topic: {title}")
            abstracts.append(f"Title: {title}\nAbstract: {abstract}")
        
        fetch_handle.close()
        return abstracts

    except Exception as e:
        print(f"  ❌ PubMed Error: {e}")
        return []
    print(f" Fetching new research from PubMed for: {query}...")
    try:
        search_handle = Entrez.esearch(db="pubmed", term=query, retmax=max_results)
        search_results = Entrez.read(search_handle)
        search_handle.close()

        id_list = search_results["IdList"]
        if not id_list: return []

        fetch_handle = Entrez.efetch(db="pubmed", id=id_list, rettype="medline", retmode="text")
        records = Medline.parse(fetch_handle)
        
        abstracts = [f"Title: {r['TI']}\nAbstract: {r['AB']}" for r in records if "AB" in r]
        fetch_handle.close()
        return abstracts
    except Exception as e:
        print(f"PubMed Error: {e}")
        return []