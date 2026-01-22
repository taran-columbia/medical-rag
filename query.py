import os
from dotenv import load_dotenv
from langchain_google_genai import GoogleGenerativeAIEmbeddings, ChatGoogleGenerativeAI
from langchain_neo4j import Neo4jGraph
from pinecone import Pinecone

load_dotenv()

# 1. Setup Gemini
embeddings = GoogleGenerativeAIEmbeddings(model="models/text-embedding-004")
llm = ChatGoogleGenerativeAI(model="gemini-flash-latest", temperature=0)

# 2. Setup Databases
graph = Neo4jGraph(
    url=os.getenv("NEO4J_URI"),
    username="neo4j",
    password=os.getenv("NEO4J_PASSWORD")
)
pc = Pinecone(api_key=os.getenv("PINECONE_API_KEY"))
index = pc.Index("medical-rag")

def medical_bot(user_question):
    print(f"\nSearching for: {user_question}")
    extraction_prompt = f"""
    Extract the 2-3 most important medical entities (drugs, diseases, symptoms) from this question.
    Return ONLY a comma-separated list of lowercase words.
    Question: {user_question}
    """

    extraction_response = llm.invoke(extraction_prompt)
    search_terms = [term.strip().lower() for term in extraction_response.content.split(",")]
    print(f"Extracted Search Terms: {search_terms}")
    
    # STEP A: Vector Search (Find relevant text)
    query_vector = embeddings.embed_query(user_question)
    vector_results = index.query(vector=query_vector, top_k=2, include_metadata=True)
    context_text = [res['metadata']['text'] for res in vector_results['matches']]
    
    # STEP B: Graph Search (Find logical connections)
    # We look for nodes mentioned in the question and find their neighbors
    graph_context = []
    for term in search_terms:
        # We use 'CONTAINS' to match partial names like 'diabetes' in 'type 2 diabetes'
        results = graph.query("""
            MATCH (n:Entity)-[r]->(m:Entity)
            WHERE n.name CONTAINS $lookup OR m.name CONTAINS $lookup
            RETURN n.name + ' ' + type(r) + ' ' + m.name AS rel
            LIMIT 3
        """, params={"lookup": term})
        graph_context.extend([r['rel'] for r in results])
    # STEP C: Combine and Answer
    combined_prompt = f"""
    You are a medical research assistant. Use the following context to answer the question.
    
    VECTOR CONTEXT (Facts): {context_text}
    GRAPH CONTEXT (Relationships): {graph_context}
    
    QUESTION: {user_question}
    ANSWER:
    """
    
    response = llm.invoke(combined_prompt)
    print("\n--- AI RESPONSE ---")
    print(response.content)


if __name__ == "__main__":
    # Test a "Multi-Hop" question!
    medical_bot("How is Metformin related to Diabetic Retinopathy?")