import streamlit as st
import scanpy as sc
import matplotlib.pyplot as plt
import pandas as pd
import io

st.title("🧬 Differential Expression (Marker Genes)")

st.markdown("""
In this step, we identify **marker genes** that are differentially expressed between clusters.  

The workflow includes:
1. **Find marker genes** – perform differential expression tests between clusters.  
   - Comparison options:  
     - **One cluster vs all other clusters**  
     - **Two specific clusters**  
   - Test method: **Wilcoxon rank-sum test** (default in Scanpy).  

2. **Inspect marker genes** – view top ranked genes and their statistics (e.g., scores, fold-changes, adjusted p-values).  

👉 Differentially expressed genes help us assign **biological meaning** to clusters.
""")


# --- Check if adata exists from Step 5 ---
if "adata" not in st.session_state:
    st.error("No AnnData object found. Please complete **Clustering & UMAP** first.")
    st.stop()

adata = st.session_state["adata"]

if "leiden" not in adata.obs:
    st.error("No clustering found. Please run **Clustering & UMAP** first.")
    st.stop()


# =========================================================
# --- Step 1: Select comparison mode ---
# =========================================================
st.subheader("📌 Step 1: Choose comparison mode")

mode = st.radio(
    "How would you like to compare clusters?",
    ["Cluster vs all other clusters", "Cluster vs cluster", "All clusters vs rest"],
    help="""
    - **Cluster vs all other clusters** 🧩: Choose one cluster and compare against all others.  
    - **Cluster vs cluster** ⚖️: Compare two specific clusters.  
    - **All clusters vs rest** 🌐: Automatically compute marker genes for *all clusters* at once.  
    """
)

clusters = sorted(adata.obs["leiden"].unique())
run_deg = False
cluster = None
cluster1 = None
cluster2 = None

if mode == "Cluster vs all other clusters":
    cluster = st.selectbox("Select cluster:", clusters)
    if st.button("Run DE analysis"):
        run_deg = True

elif mode == "Cluster vs cluster":
    col1, col2 = st.columns(2)
    with col1:
        cluster1 = st.selectbox("Cluster 1:", clusters, index=0)
    with col2:
        cluster2 = st.selectbox("Cluster 2:", clusters, index=1)
    if st.button("Run DE analysis"):
        run_deg = True

elif mode == "All clusters vs rest":
    if st.button("Run DE analysis for all clusters"):
        run_deg = True


# =========================================================
# --- Step 2: Run DE analysis ---
# =========================================================
if run_deg:
    wait_msg = st.empty()
    wait_msg.info("⏳ Running differential expression analysis...")

    if mode == "Cluster vs all other clusters":
        sc.tl.rank_genes_groups(
            adata,
            groupby="leiden",
            groups=[cluster],
            reference="rest",
            method="wilcoxon"
        )
        st.success(f"✅ Marker genes for cluster {cluster} vs all others computed.")

    elif mode == "Cluster vs cluster":
        sc.tl.rank_genes_groups(
            adata,
            groupby="leiden",
            groups=[cluster1],
            reference=cluster2,
            method="wilcoxon"
        )
        st.success(f"✅ Marker genes for cluster {cluster1} vs cluster {cluster2} computed.")

    elif mode == "All clusters vs rest":
        sc.tl.rank_genes_groups(
            adata,
            groupby="leiden",
            reference="rest",
            method="wilcoxon"
        )
        st.success("✅ Marker genes computed for ALL clusters vs rest.")

    wait_msg.empty()
    st.session_state["adata"] = adata


# =========================================================
# --- Step 2 conti.: Show results ---
# =========================================================
if "rank_genes_groups" in adata.uns:
    st.subheader("📊 Step 2: Inspect marker genes")

    st.markdown("Here are the **top ranked marker genes** per cluster (Scanpy visualization):")
    sc.pl.rank_genes_groups(adata, n_genes=20, sharey=False, show=False)
    fig = plt.gcf()
    st.pyplot(fig)
    plt.close(fig)

    # Convert results to DataFrame
    result = adata.uns["rank_genes_groups"]
    groups = result["names"].dtype.names
    dfs = []
    for g in groups:
        df = pd.DataFrame({
            "names": result["names"][g],
            "scores": result["scores"][g],
            "logfoldchanges": result["logfoldchanges"][g],
            "pvals": result["pvals"][g],
            "pvals_adj": result["pvals_adj"][g],
        })
        df["cluster"] = g
        dfs.append(df)
    df_out = pd.concat(dfs)

    st.download_button(
        label="💾 Download DE results (.csv)",
        data=df_out.to_csv(index=False).encode("utf-8"),
        file_name="DE_results.csv",
        mime="text/csv"
    )

# =========================================================
# --- Step 3: Automatically detect marker genes ---
# =========================================================
if "rank_genes_groups" in adata.uns:
    st.subheader("✨ Step 3: Automatic Marker Gene Detection")

    result = adata.uns["rank_genes_groups"]
    groups = result["names"].dtype.names

    # Collect top N marker genes per cluster
    top_n = st.number_input("Number of top genes per cluster", min_value=1, max_value=50, value=5, step=1)

    marker_dict = {}
    for g in groups:
        marker_dict[g] = result["names"][g][:top_n].tolist()

    # Convert to DataFrame for display
    marker_table = []
    for cluster, genes in marker_dict.items():
        marker_table.append({
            "Cluster": cluster,
            "Markers": ", ".join(genes)
        })
    df_markers = pd.DataFrame(marker_table)

    st.dataframe(df_markers)

    # Save marker gene list
    st.download_button(
        label="💾 Download marker genes (.csv)",
        data=df_markers.to_csv(index=False).encode("utf-8"),
        file_name="marker_genes.csv",
        mime="text/csv"
    )

    st.info("💡 These top marker genes are extracted from the DE results automatically. You can use them for cell type annotation in the next step.")


# --- Save top marker per cluster ---
if "rank_genes_groups" in adata.uns:
    result = adata.uns["rank_genes_groups"]
    groups = result["names"].dtype.names
    
    top_markers = {}
    for g in groups:
        if len(result["names"][g]) > 0:
            top_markers[g] = result["names"][g][0]  # take top 1 per cluster
    
    st.session_state["top_markers"] = top_markers
    st.info(f"💡 Saved top marker genes per cluster: {list(top_markers.values())}")


# --- Show "Next" link only after DE is done ---
    st.markdown("""
    <style>
    section[data-testid="stMain"] [data-testid="stPageLink"] a,
    section[data-testid="stMain"] [data-testid="stPageLink"] p {
      font-style: italic !important;
    }
    </style>
    """, unsafe_allow_html=True)

    spacer, right = st.columns([0.6, 0.255], gap="small")
    with right:
        st.page_link("pages/6_Assign_Cell_Type_Identity.py", label="➡️ Next: Assign Cell Identity")
