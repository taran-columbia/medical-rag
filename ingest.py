import os
import json
from dotenv import load_dotenv
from langchain_google_genai import GoogleGenerativeAIEmbeddings, ChatGoogleGenerativeAI
from langchain_neo4j import Neo4jGraph
from pinecone import Pinecone
import hashlib
import time

load_dotenv()

# 1. Setup Gemini
# embeddings = GoogleGenerativeAIEmbeddings(
#     model="models/text-embedding-004"
# )
# llm = ChatGoogleGenerativeAI(model="gemini-flash-latest", temperature=0)

# print("\n Starting ingestion process...", os.getenv("NEO4J_URI"))
# # 2. Setup Databases
# graph = Neo4jGraph(
#     url=os.getenv("NEO4J_URI"),
#     username="neo4j",
#     password=os.getenv("NEO4J_PASSWORD")
# )

# pc = Pinecone(api_key=os.getenv("PINECONE_API_KEY"))
# index = pc.Index("medical-rag")

# medical_data = [
#     "Metformin is used to treat Type 2 Diabetes by improving insulin sensitivity.",
#     "Aspirin inhibits COX-1 and COX-2 enzymes to reduce inflammation.",
#     "Type 2 Diabetes can lead to complications like Diabetic Retinopathy."
# ]

def add_in_graph(llm_response, graph):
    try:
        if isinstance(llm_response.content, str):
            raw_content = llm_response.content
        elif isinstance(llm_response.content, list):
            # Extract text from the first element of the list
            raw_content = llm_response.content[0].get('text', '')
        else:
            raw_content = str(llm_response.content)
        # Clean the response and parse JSON
        clean_json = raw_content.replace("```json", "").replace("```", "").strip()
        triples = json.loads(clean_json)

        for triple in triples:
            # We use MERGE so we don't create duplicate entities
            # We use apoc.create.relationship to handle dynamic relationship names from the LLM
            graph.query("""
                MERGE (s:Entity {name: $sub})
                MERGE (o:Entity {name: $obj})
                WITH s, o
                CALL apoc.create.relationship(s, $rel, {}, o) YIELD rel
                RETURN rel
            """, params={
                "sub": triple['sub'].strip().lower(), 
                "obj": triple['obj'].strip().lower(), 
                "rel": triple['rel'].replace(" ", "_").upper()
            })
            print(f"  Added to Graph: ({triple['sub']}) -[{triple['rel']}]-> ({triple['obj']})")
    except Exception as e:
        print(f"  Error parsing triples: {e}")

def process_data(data_list, embeddings, index, graph, llm):
    for text in data_list:
        text_id = hashlib.md5(text.encode()).hexdigest()
        print(f"Processing Chunk {text_id}: {text}")
        
        # A. Vector Step
        vector = embeddings.embed_query(text)

        # 2. DEBUG: Check the dimension before sending to Pinecone
        print(f"DEBUG: Vector length is {len(vector)}") 
        
        
        index.upsert(vectors=[{
            "id": text_id, 
            "values": vector, 
            "metadata": {"text": text}
        }])
        
        # B. Graph Step
        prompt = f"""
        Extract medical triples (Subject, Relation, Object) from this text: {text}.
        Format the output strictly as a JSON list of objects with keys 'sub', 'rel', and 'obj'.
        Do not include any other text.
        """
        response = llm.invoke(prompt)
        add_in_graph(response, graph)
        # print("  ⏳ Politeness delay (3s)...")
        # time.sleep(3)

    print("\n Successfully populated Vector and Graph databases!")

# if __name__ == "__main__":
#     process_data(medical_data)