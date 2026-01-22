import gradio as gr
from query import medical_bot
import time

def respond(question):
    # 1. Update UI to show we are starting
    yield "🔍 Initializing search and checking local database..."
    
    try:
        # We call your medical_bot logic here
        # To make it truly 'stream' status, you'd need to modify medical_bot 
        # to be a generator, but for now, let's wrap the logic:
        
        result = medical_bot(question)
        
        # 2. Return the final answer
        yield result
        
    except Exception as e:
        yield f" Error: {str(e)}"

# Define the Gradio UI
with gr.Blocks(theme=gr.themes.Soft()) as demo:
    gr.Markdown("# 🏥 MedRag: AI Medical Research Assistant")
    
    with gr.Row():
        with gr.Column():
            msg = gr.Textbox(label="Ask a Medical Question", placeholder="e.g., How does Metformin affect the heart?")
            submit = gr.Button("Run Analysis", variant="primary")
        
    with gr.Column():
        # Using Markdown here allows the AI to use bolding and lists
        output = gr.Markdown(label="Status / Analysis")

    # Link the button and the Enter key
    # show_progress="minimal" adds a subtle loading bar at the top
    submit.click(fn=respond, inputs=msg, outputs=output, show_progress="full")
    msg.submit(fn=respond, inputs=msg, outputs=output, show_progress="full")

if __name__ == "__main__":
    demo.launch(share=True)