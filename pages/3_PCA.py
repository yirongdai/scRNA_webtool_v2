import streamlit as st
import scanpy as sc
import matplotlib.pyplot as plt
import math

st.title("📉 Linear Dimensional Reduction (PCA)")

st.markdown("""
Principal Component Analysis (**PCA**) is the first dimensional reduction step in the standard scRNA-seq workflow.  
It projects cells into a low-dimensional space that captures the **major sources of variation** in the data.

The workflow includes:  
1. **Run PCA** – choose how many principal components to compute and whether to restrict the analysis to **highly variable genes (HVGs)**.  
2. **Inspect variance explained (Elbow plot)** – determine how many PCs capture most of the variation.  
3. **Visualize PCA scatter plots** – explore how cells separate in the PCA space.  
4. **Examine PCA loadings** – identify which genes contribute most to each principal component.  

👉 *These PCs will be used in downstream steps such as clustering and UMAP.*
""")

# --- Check if adata exists ---
if "adata" not in st.session_state:
    st.error("No AnnData object found. Please complete **Preprocessing** first.")
    st.stop()

adata = st.session_state["adata"]

# --- User options ---
n_comps = st.number_input(
    "Number of PCs", 
    min_value=5, max_value=100, value=50, step=5,
    help="How many principal components to compute. Typically 30–50."
)

use_hvg_option = st.selectbox(
    "Use highly variable genes (HVGs)?", 
    options=["Auto (use HVGs if available)", "Yes (only HVGs)", "No (all genes)"],
    index=0,
    help="Auto – If HVGs have been identified, PCA will use only HVGs. "
        "If not, PCA will use all genes.\n\n"
        "Yes – Force PCA to use only HVGs.\n"
        "No – Force PCA to use all genes."
)

if use_hvg_option == "Auto (default)":
    use_highly_variable = None
elif use_hvg_option == "Yes (only HVGs)":
    use_highly_variable = True
else:
    use_highly_variable = False

# --- Run PCA ---
if st.button("Run PCA"):
    wait_placeholder = st.empty()
    wait_placeholder.info("⏳ Running PCA...")

    sc.tl.pca(
        adata,
        n_comps=n_comps,
        use_highly_variable=use_highly_variable,
        svd_solver="arpack"
    )

    # ✅ Save selection immediately
    adata.uns["n_pcs_selected"] = n_comps
    st.session_state["n_pcs"] = n_comps

    st.success(f"✅ PCA done with n_comps={n_comps}")
    st.session_state["adata"] = adata
    wait_placeholder.empty()


# --- If PCA done, show plots ---
if "X_pca" in adata.obsm_keys():
    # ---- Elbow plot ----
    st.subheader("📊 Variance explained (Elbow plot)")
    st.markdown("Use this plot to estimate how many PCs capture most of the variance.")

    n_pcs_elbow = st.number_input(
        "Number of PCs to display:",
        min_value=5, max_value=100, value=20, step=5
    )
    sc.pl.pca_variance_ratio(adata, log=False, n_pcs=n_pcs_elbow, show=False)
    fig = plt.gcf()
    st.pyplot(fig)
    plt.close(fig)

    # --- User chooses PCs for downstream analysis ---
    st.subheader("🔧 Select number of PCs for downstream analysis")
    st.markdown("""This is based on where the elbow plot plateaus.""")
    n_pcs_final = st.number_input(
        "Number of PCs to use (for clustering, UMAP, etc.):",
        min_value=5, max_value=100, value=30, step=5
    )

    if st.button("Save selection"):
        st.session_state["n_pcs"] = n_pcs_final
        st.success(f"✅ Using {n_pcs_final} PCs for downstream analysis.")

        # Save AnnData with PCA + chosen n_pcs stored in uns
        adata.uns["n_pcs_selected"] = n_pcs_final
        import tempfile, os
        with tempfile.NamedTemporaryFile(delete=False, suffix=".h5ad") as tmp:
            adata.write(tmp.name)
            tmp_path = tmp.name
        with open(tmp_path, "rb") as f:
            st.download_button(
                "💾 Download PCA-processed data (.h5ad)",
                f,
                "adata_pca.h5ad",
                "application/octet-stream"
            )
        os.remove(tmp_path)

    # ---- PCA scatter plot ----
    st.subheader("🖼️ PCA Scatter Plot")
    st.markdown("""
                Visualize cells in the PCA space. Choose PCs for X and Y axes.  
                You can color cells by **gene expression** to see which clusters are associated with certain marker genes.
                Each point represents a **single cell**.  
                - **Axes (PC1, PC2, …)**: principal components capturing major sources of variation.  
                - **Colors**: expression level of the selected gene in each cell  
                    - Dark = low or no expression  
                    - Bright = high expression  
                """)
   
    
    st.info("""
    💡 **Tip:** Marker genes are genes whose expression highlights specific cell types.  
    Here are some commonly used marker genes in PBMC data:

    - **CST3** → dendritic cell / monocyte marker  
    - **NKG7** → NK cell / cytotoxic T cell marker  
    - **MS4A1** → B cell marker  
    - **CD3D** → T cells  
    - **PPBP** → Platelets  
    - **S100A4** → CD4 memory T cells
    """)

    # Choosing marker genes
    marker_genes = ["CST3", "NKG7", "MS4A1", "CD3D", "PPBP", "S100A4"]

    # Input number of PCs
    x_pc = st.number_input("PC for X-axis", min_value=1, max_value=50, value=1, step=1)
    y_pc = st.number_input("PC for Y-axis", min_value=1, max_value=50, value=2, step=1)

    # Select genes
    selected_genes = st.multiselect(
        "Color cells by gene(s) (optional):",
        options=adata.var_names.tolist(),
        default=marker_genes,
        help="Choose one or more genes. If left empty, no gene coloring will be applied."
    )

    # Plotting PCA
    if st.button("Plot PCA Scatter"):
        sc.pl.pca(
            adata,
            components=f"{x_pc},{y_pc}",
            color=selected_genes if selected_genes else None,
            show=False
        )
        fig = plt.gcf()
        st.pyplot(fig)
        plt.close(fig)


    # ---- PCA Loadings ----
    st.subheader("🔎 Top contributing genes per PC")
    st.markdown("See which genes drive each principal component (PC).")

    n_pcs_loadings = st.number_input(
        "Number of PCs to inspect (must be even number to display):",
        min_value=2, max_value=20, value=6, step=2
    )

    n_pairs = math.ceil(n_pcs_loadings / 2)
    for i in range(n_pairs):
        pc1 = 2 * i + 1
        pc2 = 2 * i + 2
        if pc2 <= n_pcs_loadings:
            sc.pl.pca_loadings(
                adata,
                components=(pc1, pc2),
                include_lowest=True,
                show=False
            )
            fig = plt.gcf()
            st.pyplot(fig)
            plt.close(fig)

else:
    st.info("👉 Run PCA first to view elbow plot, scatter plot, and PC loadings.")

# --- Show "Next: Clustering & UMAP" only after PCA is done ---
if "X_pca" in adata.obsm_keys():
    # italicize page_link labels in the main area
    st.markdown("""
    <style>
    section[data-testid="stMain"] [data-testid="stPageLink"] a,
    section[data-testid="stMain"] [data-testid="stPageLink"] p {
      font-style: italic !important;
    }
    </style>
    """, unsafe_allow_html=True)

    # push the link to the right; increase left column width to push link further right
    spacer, right = st.columns([0.5, 0.2], gap="small")
    with right:
        st.page_link("pages/4_Clustering & UMAP.py", label="➡️ Next: Clustering & UMAP")

