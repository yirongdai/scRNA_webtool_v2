# scRNA-seq Analysis Webtool

This project is part of my **PhD coursework**, where I am developing an interactive [Streamlit](https://streamlit.io/) app for processing **single-cell RNA-seq (scRNA-seq) data**.

The tool is inspired by workflows from:  
- [Seurat PBMC3k tutorial](https://satijalab.org/seurat/articles/pbmc3k_tutorial) from Satija Lab  
- [Scanpy PBMC3k tutorial](https://scanpy-tutorials.readthedocs.io/en/latest/pbmc3k.html#preprocessing)  

It allows users to perform essential scRNA-seq analysis **step by step** through a simple web interface.

> This project was created with the help of **ChatGPT-5** (OpenAI), for guidance on code structure, deployment, and documentation.

---

## Try it out
[Launch the webtool here](https://yirongdai-scrna-webtool-home-j2ukrh.streamlit.app/)  

---

## Data Input Options

You can start analysis in three ways:

1. **Upload a preprocessed file**  
   - `.h5ad` (AnnData format)  
   - Recommended for large datasets  

2. **Upload raw 10X files (individually)**  
   - `matrix.mtx` or `matrix.mtx.gz`  
   - `genes.tsv` or `features.tsv.gz`  
   - `barcodes.tsv` or `barcodes.tsv.gz`  
   - These will be converted internally into `.h5ad`  

3. **Use Demo Data**  
   - PBMC3k dataset from 10X Genomics  
   - Download link: [pbmc3k_filtered_gene_bc_matrices.tar.gz](https://cf.10xgenomics.com/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz)  
   - Already included in the `data/` folder for quick testing

---

## Workflow

The app guides you through the following steps:

1. **Load Data**  
2. **QC Filtering** – remove low-quality cells, high mito %  
3. **Preprocessing** – normalization & highly variable gene selection  
4. **Run PCA** – dimensionality reduction  
5. **Clustering** – Leiden/Louvain algorithms  
6. **Run UMAP** – visualize cells in 2D  
7. **DEGs** – marker gene identification  
8. **Assign Cell Type Identity** – annotation with known markers  

---

## Installation (Local Development)

Clone the repo and install dependencies:

```bash
git clone https://github.com/yirongdai/scRNA_webtool.git
cd scRNA_webtool
pip install -r requirements.txt
