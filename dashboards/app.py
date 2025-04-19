import streamlit as st, pandas as pd, plotly.express as px
from sklearn.manifold import TSNE
from agents import analyst
from config import settings

ABSTRACTS_STORE, PDF_STORE = "abstracts_store_id", "pdf_store_id"

st.title("Literag – RAG Dashboard")

tab1, tab2 = st.tabs(["Search", "Embedding map"])

with tab1:
    q = st.text_input("Ask a question")
    if st.button("Submit") and q:
        res = analyst.search(q, PDF_STORE)
        st.markdown("### Top chunks")
        for i, r in enumerate(res, 1):
            st.markdown(f"**{i}.** {r['text']}  \nDOI: `{r['metadata']['doi']}`")

with tab2:
    st.markdown("UMAP of abstract embeddings")
    # In practice load embeddings from store → X
    X = pd.read_pickle("abstract_embeds.pkl")
    coords = TSNE(n_components=2).fit_transform(X)
    fig = px.scatter(x=coords[:,0], y=coords[:,1])
    st.plotly_chart(fig, use_container_width=True)
