import streamlit as st

st.set_page_config(page_title="scRNA-seq Analysis Webtool", page_icon="🧬")

st.title("🧬 scRNA-seq Analysis Webtool")

st.markdown("""
Welcome to the **Single-cell RNA-seq Analysis Webtool** 👋  
This interactive app was developed as part of my **PhD coursework** to make
the [Scanpy](https://scanpy.readthedocs.io/) and [Seurat](https://satijalab.org/seurat/) 
workflows more accessible through a simple, step-by-step web interface.

The tool is inspired by:
- [Seurat PBMC3k tutorial](https://satijalab.org/seurat/articles/pbmc3k_tutorial) (Satija Lab)  
- [Scanpy PBMC3k tutorial](https://scanpy-tutorials.readthedocs.io/en/latest/pbmc3k.html)  

---

## 🚀 Workflow Overview

You can explore your single-cell dataset using the following steps:

1. **Load Data** – Upload `.h5ad`, `.h5`, or `.loom` files, raw 10X files, or use demo PBMC3k data  
2. **Preprocessing** – Perform QC filtering, normalization, HVG selection, and scaling  
3. **PCA** – Linear dimensional reduction to capture major sources of variation  
4. **Clustering & UMAP** – Construct neighborhood graph, cluster cells (Leiden/Louvain), and embed with UMAP  
5. **DEGs** – Find marker genes with differential expression tests  
   - One cluster vs all other clusters  
   - Two specific clusters  
6. **Assign Cell Type Identity** – Explore marker expression and auto-annotate with **CellTypist**  

---

## 📦 Data Input Options
- Upload `.h5ad`, `.h5`, or `.loom` file (recommended for large datasets)  
- Upload **raw 10X files** (`matrix.mtx`, `genes.tsv/features.tsv`, `barcodes.tsv`) → auto-converts to `.h5ad`  
- Use the included **PBMC3k demo dataset**  

---

ℹ️ This project was built with the help of **ChatGPT-5 (OpenAI)** for code structure, deployment, and documentation.
""")
