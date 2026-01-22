import os
from dotenv import load_dotenv
from langchain_google_genai import GoogleGenerativeAIEmbeddings, ChatGoogleGenerativeAI
from langchain_neo4j import Neo4jGraph
from pinecone import Pinecone
from fetcher import fetch_medical_papers
from ingest import process_data

load_dotenv()

# 1. Setup Gemini
embeddings = GoogleGenerativeAIEmbeddings(model="models/text-embedding-004")
llm = ChatGoogleGenerativeAI(model="gemini-flash-latest", temperature=0).with_retry(stop_after_attempt=6)

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
    vector_results = index.query(vector=query_vector, top_k=1, include_metadata=True)
    top_score = vector_results['matches'][0]['score'] if vector_results['matches'] else 0

    if top_score < 0.65:
        print(f"Low confidence ({top_score:.2f}). Accessing PubMed...")
        new_data = fetch_medical_papers(extraction_response.content, max_results=3)
        print(f"Fetched {len(new_data)} new papers.")
        if new_data:
            # Pass our initialized objects to ensure consistency
            process_data(new_data, embeddings, index, graph, llm)
            # Re-run search
            vector_results = index.query(vector=query_vector, top_k=2, include_metadata=True)
    # STEP B: Graph Search (Find logical connections)
    # We look for nodes mentioned in the question and find their neighbors

    
    graph_context = []
    for term in search_terms:
        # We use 'CONTAINS' to match partial names like 'diabetes' in 'type 2 diabetes'
        cypher_results = graph.query("""
            MATCH (n:Entity)-[r]->(m:Entity)
            WHERE n.name CONTAINS $lookup OR m.name CONTAINS $lookup
            RETURN n.name + ' ' + type(r) + ' ' + m.name AS rel
            LIMIT 3
        """, params={"lookup": term})
        graph_context.extend([r['rel'] for r in cypher_results])
    # STEP C: Combine and Answer
    context_text = [res['metadata']['text'] for res in vector_results['matches']]
    combined_prompt = f"""
    You are a professional medical research assistant. 

    1. If the provided VECTOR and GRAPH CONTEXT contain the answer, prioritize that information.
    2. If the context is missing specific details (like common causes), you MAY use your general medical knowledge to provide a complete answer.
    3. ALWAYS distinguish between what was found in the research papers vs. general medical facts.

    CONTEXT:
    Text Facts: {context_text}
    Relationships: {graph_context}

    QUESTION: {user_question}
    """
    
    response = llm.invoke(combined_prompt)
    if isinstance(response.content, str):
        final_answer = response.content
    elif isinstance(response.content, list):
        # This handles the specific list format you saw earlier
        final_answer = response.content[0].get('text', str(response.content))
    else:
        final_answer = str(response.content)

    print("\n--- AI RESPONSE ---")
    print(final_answer)
    return final_answer
    


if __name__ == "__main__":
    # Test a "Multi-Hop" question!
    medical_bot("What is some common reasons for hypertension?")