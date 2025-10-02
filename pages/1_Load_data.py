import streamlit as st
import scanpy as sc
import tempfile, os, io
import pandas as pd
import numpy as np

# --- Sidebar controls ---
st.sidebar.header("Step 1: Load Data 📂")

option = st.sidebar.radio(
    "Choose how to load your data:",
    ["Upload files (.h5ad / .h5 / .loom)",
     "Upload 10X files (matrix + genes/features + barcodes)", 
     "Use Demo Data"],
    index=2
)

adata = None

# 1. Upload h5ad, h5, or loom
if option == "Upload files (.h5ad / .h5 / .loom)":
    uploaded_files = st.sidebar.file_uploader(
        "Upload one or more files",
        type=["h5ad", "h5", "loom"],
        accept_multiple_files=True
    )

    if uploaded_files:
        adata_list = []
        for uploaded_file in uploaded_files:
            try:
                with tempfile.TemporaryDirectory() as tmpdir:
                    tmp_path = os.path.join(tmpdir, uploaded_file.name)
                    with open(tmp_path, "wb") as f:
                        f.write(uploaded_file.read())

                    # Load depending on suffix
                    if uploaded_file.name.endswith(".h5ad"):
                        adata_temp = sc.read_h5ad(tmp_path)
                    elif uploaded_file.name.endswith(".loom"):
                        adata_temp = sc.read_loom(tmp_path)
                    elif uploaded_file.name.endswith(".h5"):
                        adata_temp = sc.read_10x_h5(tmp_path)
                    else:
                        st.sidebar.error(f"❌ Unsupported file type: {uploaded_file.name}")
                        continue

                    adata_list.append(adata_temp)
                    st.sidebar.success(f"✅ Loaded {uploaded_file.name}")

            except Exception as e:
                st.sidebar.error(f"❌ Error reading {uploaded_file.name}: {e}")

        # Merge logic
        if len(adata_list) == 1:
            adata = adata_list[0]
        elif len(adata_list) > 1:
            try:
                adata = adata_list[0].concatenate(*adata_list[1:], batch_key="batch")
                st.sidebar.success("✅ Multiple files merged successfully with `batch` column.")
            except Exception as e:
                st.sidebar.error(f"❌ Failed to merge files: {e}")


# 2. Upload 10X files individually
elif option == "Upload 10X files (matrix + genes/features + barcodes)":
    uploaded_files = st.sidebar.file_uploader(
        "Upload your 10X files (matrix.mtx + genes.tsv/features.tsv + barcodes.tsv)",
        type=["mtx", "tsv", "gz"],
        accept_multiple_files=True
    )

    if uploaded_files:
        # classify by name
        mtx_file = next((f for f in uploaded_files if f.name.endswith(".mtx") or f.name.endswith(".mtx.gz")), None)
        genes_file = next((f for f in uploaded_files if "genes" in f.name or "features" in f.name), None)
        barcodes_file = next((f for f in uploaded_files if "barcodes" in f.name), None)

        if mtx_file and genes_file and barcodes_file:
            with tempfile.TemporaryDirectory() as tmpdir:
                try:
                    # Decide target filenames
                    mtx_name = "matrix.mtx.gz" if mtx_file.name.endswith(".gz") else "matrix.mtx"
                    if genes_file.name.endswith(".gz"):
                        feat_name = "features.tsv.gz" if "features" in genes_file.name else "genes.tsv.gz"
                    else:
                        feat_name = "features.tsv" if "features" in genes_file.name else "genes.tsv"
                    bar_name = "barcodes.tsv.gz" if barcodes_file.name.endswith(".gz") else "barcodes.tsv"

                    # Save files
                    for f, fname in [(mtx_file, mtx_name), (genes_file, feat_name), (barcodes_file, bar_name)]:
                        with open(os.path.join(tmpdir, fname), "wb") as out_f:
                            out_f.write(f.read())

                    # Read with Scanpy
                    try:
                        adata = sc.read_10x_mtx(tmpdir, var_names="gene_symbols", cache=False)
                    except Exception:
                        adata = sc.read_10x_mtx(tmpdir, var_names="gene_ids", cache=False)

                    st.sidebar.success("✅ Data loaded from 10X files!")

                except Exception as e:
                    st.sidebar.error(f"❌ Error reading 10X files: {e}")


# 3. Use demo data
elif option == "Use Demo Data":
    try:
        adata = sc.datasets.pbmc3k()
        st.sidebar.success("✅ Demo data loaded (PBMC 3k)!")
    except Exception as e:
        st.sidebar.error(f"❌ Could not load demo dataset: {e}")


# --- Main page ---
st.title("📂 Load your data here!")

st.markdown("""
Welcome to the **data loading step** of the scRNA-seq Webtool!  
Choose one of the three options in the sidebar to load your dataset.
""")

# Store and preview
if adata is not None:
    # Make gene names unique to avoid duplicate column errors
    adata.var_names_make_unique()

    # Save to session state
    st.session_state["adata"] = adata

    st.subheader("🔍 Preview")
    st.write(f"**Number of cells:** {adata.n_obs}")
    st.write(f"**Number of genes:** {adata.n_vars}")

    try:
        X = adata.X[:5, :5].toarray() if hasattr(adata.X, "toarray") else adata.X[:5, :5]
        df = pd.DataFrame(X, index=adata.obs_names[:5], columns=adata.var_names[:5])
        st.dataframe(df)
    except Exception as e:
        st.error(f"Could not display preview matrix: {e}")
else:
    st.info("Please upload a dataset or select demo data to continue.")

# --- Bottom-right "Next: Preprocessing" link ---
# Right-align & italicize the "Next: Preprocessing" page link in the main area
st.markdown("""
<style>
/* Only affect links rendered in the main content, not the sidebar */
section[data-testid="stMain"] [data-testid="stPageLink"] a,
section[data-testid="stMain"] [data-testid="stPageLink"] p {
  font-style: italic !important;
}
</style>
""", unsafe_allow_html=True)


spacer, right = st.columns([0.6, 0.2], gap="small")


# 2) Visibility logic:
#    - For Option 3 (demo): show by default
#    - For Option 1/2: show only when data has been uploaded (adata is not None)
show_next = (option == "Use Demo Data") or (adata is not None)
if show_next:
    with right:
        st.page_link(
            "pages/2_Preprocessing.py",
            label="➡️ Next: Preprocessing",
        )