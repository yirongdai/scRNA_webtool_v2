import streamlit as st
import scanpy as sc
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import io


st.title("🧭 Gene Expression and Cell Type Annotation")

st.markdown("""
In this step, we explore **gene expression patterns** on UMAP and annotate clusters with **cell type labels**.

The workflow includes:
1. **Gene expression on UMAP** – visualize expression of marker or user-selected genes.  
2. **Automatic marker gene detection** – compute top genes distinguishing clusters.  
3. **Cell type annotation (CellTypist)** – automatically predict cell types using pretrained models.  

👉 This step helps us assign **biological meaning** to the clusters identified earlier.
""")

# --- Check if adata exists from Step 6 ---
if "adata" not in st.session_state:
    st.error("No AnnData object found. Please complete **DEGs** first.")
    st.stop()

adata = st.session_state["adata"]

# =========================================================
# Part 1: Gene Expression on UMAP
# =========================================================
st.subheader("📌 Step 1: Gene Expression on UMAP")

# Marker gene tip
st.info("""
💡 **Tip:** Marker genes are genes whose expression highlights specific cell types.  
Here are some commonly used marker genes in PBMC data:

- **CST3** → dendritic cell / monocyte marker  
- **NKG7** → NK cell / cytotoxic T cell marker  
- **MS4A1** → B cell marker  
- **CD3D** → T cell marker  
- **PPBP** → Platelet marker 
- **S100A4** → CD4 memory T cell marker
""")

# Default marker list
default_markers = ["CST3", "NKG7", "MS4A1", "CD3D", "PPBP", "S100A4"]

# Gene list source
gene_source = st.radio(
    "Choose gene list:",
    ["Highly variable genes", "All genes"],
    index=0,
    help="""
- **Highly variable genes (HVGs)**: Focus only on the most informative genes (faster, more concise list).  
- **All genes**: Full list of genes in the dataset. Choose this if your marker gene is not in HVGs.  
"""
)

# Build gene list
all_genes = list(adata.var_names)
if "highly_variable" in adata.var.columns:
    hvg_genes = list(adata.var[adata.var["highly_variable"]].index)
else:
    hvg_genes = []

gene_list = hvg_genes if (gene_source == "Highly variable genes" and len(hvg_genes) > 0) else all_genes

# Keep only defaults that exist in gene_list
valid_defaults = [g for g in default_markers if g in gene_list]

# Select genes
selected_genes = st.multiselect(
    "Select one or more genes to visualize:",
    options=gene_list,
    default=valid_defaults,
    help="Choose from the list of genes to color cells on UMAP."
)

st.markdown("""
You can explore how genes are expressed across clusters using different visualization methods:  
- **UMAP** → Leiden clusters + selected gene(s) expression (side by side).  
- **Violin plots** → Distribution of expression across clusters.  
""")

plot_type = st.radio(
    "Choose visualization type:",
    ["UMAP", "Violin plots"]
)

# --- UMAP option ---
if plot_type == "UMAP":
    if st.button("Plot UMAP with selected genes"):
        st.session_state["last_selected_genes"] = selected_genes
        st.session_state["show_gene_umaps"] = True

    if st.session_state.get("show_gene_umaps", False) and st.session_state.get("last_selected_genes"):
        colors = ["leiden"] + st.session_state["last_selected_genes"]

        # Show plots 2 per row
        for i in range(0, len(colors), 2):
            cols = st.columns(2)
            for j in range(2):
                if i + j < len(colors):
                    gene = colors[i + j]
                    with cols[j]:
                        st.subheader(f"UMAP: {gene}")
                        sc.pl.umap(adata, color=gene, show=False, use_raw=False)
                        fig = plt.gcf()
                        st.pyplot(fig)
                        plt.close(fig)

# --- Violin plot option ---
elif plot_type == "Violin plots":
    if st.button("Plot violin plots"):
        for gene in selected_genes:
            st.subheader(f"Violin plot: {gene}")
            sc.pl.violin(adata, keys=gene, groupby="leiden", show=False)
            fig = plt.gcf()
            st.pyplot(fig)
            plt.close(fig)

# =========================================================
# Part 2: Automatic Marker Gene Detection
# =========================================================
st.subheader("📌 Step 2: Automatic Marker Gene Detection")

st.markdown("""
Here we use **differential expression analysis** to automatically find **marker genes** for each cluster.  
These are genes that are **highly expressed in one cluster compared to others**.
""")

top_n = st.slider("Number of top marker genes per cluster:", min_value=3, max_value=20, value=5)

if st.button("Find marker genes"):
    sc.tl.rank_genes_groups(
        adata,
        groupby="leiden",
        method="wilcoxon"
    )
    st.success(f"✅ Computed marker genes for all clusters (top {top_n}).")
    st.session_state["adata"] = adata

    # Plot Scanpy result
    sc.pl.rank_genes_groups(adata, n_genes=top_n, sharey=False, show=False)
    fig = plt.gcf()
    st.pyplot(fig)
    plt.close(fig)

    # Convert to DataFrame
    result = adata.uns["rank_genes_groups"]
    groups = result["names"].dtype.names
    dfs = []
    for g in groups:
        df = pd.DataFrame({
            "names": result["names"][g][:top_n],
            "scores": result["scores"][g][:top_n],
            "logfoldchanges": result["logfoldchanges"][g][:top_n],
            "pvals_adj": [f"{x:.2e}" for x in result["pvals_adj"][g][:top_n]]
        })
        df["cluster"] = g
        dfs.append(df)
    df_out = pd.concat(dfs)

    st.dataframe(df_out)

    st.download_button(
        label="💾 Download marker genes (.csv)",
        data=df_out.to_csv(index=False).encode("utf-8"),
        file_name="marker_genes.csv",
        mime="text/csv"
    )

# =========================================================
# Part 3: Cell Type Annotation (CellTypist)
# =========================================================
st.subheader("📌 Step 3: Cell Type Annotation (CellTypist)")

st.markdown("""
So far, clustering has only given us **cluster numbers** (0, 1, 2 …).  
These numbers show groups of similar cells, but they don’t tell us **what type of cells** they are.

Here we use **CellTypist**, a machine learning tool trained on thousands of annotated single-cell datasets,  
to automatically predict the **biological identity** of each cell (e.g., T cells, B cells, NK cells, monocytes, platelets).

👉 In short:  
- **Leiden clustering** = mathematical grouping of similar cells.  
- **CellTypist** = translate those groups into known **cell types**.  

The results will be shown on the UMAP plot and as a summary table.
""")



try:
    import celltypist
    import pandas as pd

    model_choice = st.selectbox(
        "Choose CellTypist model:",
        ["Immune_All_High.pkl", "Immune_All_Low.pkl"]
    )

    if st.button("Run CellTypist annotation"):
        wait_placeholder = st.empty()
        wait_placeholder.info("⏳ Please wait a moment while running CellTypist...")

        # Load model
        model = celltypist.models.Model.load(model_choice)

        # Run on raw normalized data (stored in Step 3)
        prediction = celltypist.annotate(
            adata.raw.to_adata(),
            model=model,
            majority_voting=True
        )

        # --- Handle CellTypist outputs safely ---
        if hasattr(prediction, "predicted_labels"):
            labels_df = prediction.predicted_labels
            if isinstance(labels_df, pd.DataFrame):
                if "majority_voting" in labels_df.columns:
                    adata.obs["predicted_labels"] = labels_df["majority_voting"]
                else:
                    adata.obs["predicted_labels"] = labels_df.iloc[:, 0]
            else:
                adata.obs["predicted_labels"] = labels_df
        elif "predicted_labels" in prediction.adata.obs:
            adata.obs["predicted_labels"] = prediction.adata.obs["predicted_labels"]
        elif "majority_voting" in prediction.adata.obs:
            adata.obs["predicted_labels"] = prediction.adata.obs["majority_voting"]
        else:
            st.error("Could not find predicted labels in CellTypist output.")
            st.stop()

        st.success(f"✅ CellTypist annotation complete using {model_choice}")
        st.session_state["adata"] = adata

        # --- UMAP with annotation ---
        st.subheader("UMAP with predicted cell types")
        sc.pl.umap(
            adata,
            color="predicted_labels",
            legend_loc="on data",
            show=False,
            use_raw=False
        )
        fig = plt.gcf()
        st.pyplot(fig)
        plt.close(fig)

        # --- Table of predicted cell types ---
        st.subheader("Cell type counts")
        st.dataframe(adata.obs["predicted_labels"].value_counts())

        # --- Download CSV of predicted labels ---
        st.subheader("Download annotations")
        df_labels = adata.obs[["predicted_labels"]].copy()
        csv_data = df_labels.to_csv().encode("utf-8")
        st.download_button(
            label="💾 Download predicted labels (.csv)",
            data=csv_data,
            file_name="celltypist_predicted_labels.csv",
            mime="text/csv"
        )

        wait_placeholder.empty()

except ImportError:
    st.error("CellTypist is not installed. Please run: `pip install celltypist`")

# =========================================================
# 📊 Step 4: Visualize marker expression by cell types
# =========================================================
st.subheader("📌 Step 4: Visualize marker expression by cell types")

# --- Check dependency ---
if "adata" not in st.session_state or "predicted_labels" not in st.session_state["adata"].obs:
    st.error("❌ Please complete Step 3 (CellTypist annotation) before visualizing.")
    st.stop()

adata = st.session_state["adata"]

# --- Marker gene source ---
marker_source = st.radio(
    "Choose marker gene source:",
    ["Custom selection", "Top 1 per cluster from DEGs"],
    index=0
)

if marker_source == "Custom selection":
    # --- Marker gene suggestions ---
    st.info("""
    💡 **Tip:** Try common PBMC marker genes:  
    - **CST3** → monocytes/dendritic cells  
    - **NKG7** → NK / cytotoxic T cells  
    - **MS4A1** → B cells  
    - **CD3D** → T cells  
    - **PPBP** → Platelets  
    - **S100A4** → CD4 memory T cells
    """)

    # Default markers
    default_markers = ["CST3", "NKG7", "MS4A1", "CD3D", "PPBP", "S100A4"]

    marker_genes = st.multiselect(
        "Select marker genes to plot:",
        options=adata.var_names.tolist(),
        default=[g for g in default_markers if g in adata.var_names]
    )
else:
    if "top_markers" in st.session_state:
        marker_genes = list(st.session_state["top_markers"].values())
        st.success(f"✅ Using top marker genes from DEGs: {marker_genes}")
        st.info("ℹ️ If you want to use **top 1 gene from all clusters**, please go to **DEGs** → "
                "**Step 1: Choose comparison mode** and select **All clusters vs rest**, then run DE analysis.")
    else:
        st.error("❌ No top markers found. Please complete DEGs first.")
        marker_genes = []


st.markdown("""
Now that we have **cell type annotations**, we can visualize marker genes across cell types.  
Choose a visualization style below:
- **Dotplot** → shows fraction of cells (dot size) and mean expression (color).  
- **Stacked violin plot** → compact view of expression distributions across groups.  
""")

plot_type = st.radio(
    "Choose plot type:",
    ["Dotplot", "Stacked violin"],
    index=0
)

if st.button("Generate plot"):
    if len(marker_genes) == 0:
        st.warning("⚠️ Please select at least one marker gene to plot.")
    else:
        buf = io.BytesIO()
        if plot_type == "Dotplot":
            sc.pl.dotplot(adata, marker_genes, groupby="predicted_labels", show=False)
        else:
            sc.pl.stacked_violin(adata, marker_genes, groupby="predicted_labels", show=False)
        
        fig = plt.gcf()
        st.pyplot(fig)


        # --- Save plot ---
        fig.savefig(buf, format="png", dpi=300, bbox_inches="tight")
        st.download_button(
            label="💾 Download plot (.png)",
            data=buf.getvalue(),
            file_name=f"{plot_type}_marker_expression.png",
            mime="image/png"
        )
        buf.close()
        plt.close(fig)
