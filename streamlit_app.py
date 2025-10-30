import streamlit as st
import pandas as pd
from utils import file_handlers
import os
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import io
from io import StringIO
from Bio import Phylo
from genome_analyzer import all_vs_all_fastani, neighbor_joining_tree
from version import __version__

st.title("GenomeCheck ✅")

# Short description (subheader-sized) directly under the title
st.subheader("Quickly assess, compare, and verify genome assemblies - all without the command line!")

# Initialize session-specific cache to prevent cross-talk between users
if 'session_id' not in st.session_state:
    import uuid
    st.session_state.session_id = str(uuid.uuid4())
    
    # Clear any existing cached data from previous sessions
    keys_to_clear = [key for key in st.session_state.keys() 
                    if key.startswith(('fastani_', 'heatmap_', 'tree_', 'ani_', 'genome_', 'qc_', 'files_processed', 'temp_files_created'))]
    for key in keys_to_clear:
        del st.session_state[key]
    # Create a session-specific temporary directory and record it in session state
    import os
    session_temp = os.path.join('./temp_sessions', st.session_state.session_id)
    st.session_state.session_temp = session_temp
    # Export to environment so library code can pick it up
    os.environ['SESSION_TEMP_DIR'] = session_temp
    # Ensure the directory exists
    file_handlers.create_temp_directory(session_temp)

# Link to repository / documentation
st.markdown(
    "<div style='background:#fff8e1;border-left:6px solid #ffd54f;padding:12px;border-radius:6px'>"
    "<h3 style='margin:0'>📘 Repository & Documentation</h2>"
    "<p style='margin:6px 0 0'>Visit <a href='https://github.com/meadm/GenomeCheck' target='_blank' rel='noopener'>the project on GitHub</a> for full docs, releases, and example data.</p>"
    "</div>",
    unsafe_allow_html=True,
)

# Version info
st.caption(f"Version: {__version__}")

# File Upload
# Small spacer before the uploader for visual separation
# Slightly smaller-than-subheader upload label rendered via Markdown
st.write("---")
st.subheader("Upload Genome FASTA Files")

# file_uploader requires a non-empty label for accessibility; hide the label visually
uploaded_files = st.file_uploader(
    "Upload genome FASTA files (hidden label)",
    type=["fasta", "fa", "fna"],
    accept_multiple_files=True,
    label_visibility="collapsed",
)

if uploaded_files:
    st.write(f"{len(uploaded_files)} genomes uploaded")
    
    # Clear previous session data to prevent cross-talk between users
    cache_keys_to_clear = [key for key in st.session_state.keys() 
                          if key.startswith(('fastani_', 'heatmap_', 'tree_', 'ani_', 'genome_', 'files_processed', 'temp_files_created'))]
    for key in cache_keys_to_clear:
        del st.session_state[key]

    # BUSCO options
    st.write("---")
    st.subheader("Genome Statistics and BUSCO Options")
    # Detect BUSCO availability in the environment
    import shutil
    busco_available = shutil.which("busco") is not None

    include_busco = st.checkbox(
        "Include BUSCO analysis (genome completeness assessment that can take several minutes per file)",
        value=False,
        help="BUSCO assesses genome completeness by searching for universal single-copy orthologs. This may take several minutes per genome.",
        disabled=not busco_available,
    )
    # If BUSCO is not available (e.g., Streamlit Cloud), force-disable and inform the user
    if not busco_available:
        st.warning("BUSCO is not available on the Streamlit Cloud or lean Docker (e.g. 'genomecheck:latest') implementations of this app. Run locally or use our Docker image with BUSCO ('genomecheck:busco') if you need BUSCO analyses.")
        include_busco = False
    
    if include_busco:
        st.info("⚠️ **Note**: BUSCO analysis can take 5-30 minutes per genome depending on size and complexity.")
        
        # Lineage selection
        busco_lineage = st.selectbox(
            "Select organism lineage for BUSCO analysis:",
            options=[
                ("bacteria_odb12", "Bacteria (odb12)"),
                ("eukaryota_odb12", "Eukaryota (odb12)"),
                ("archaea_odb12", "Archaea (odb12)"),
                ("custom", "Custom lineage (enter manually)")
            ],
            format_func=lambda x: x[1],
            help="Choose the appropriate lineage for your organisms."
        )
        
        # Extract the lineage code or allow custom input
        if busco_lineage[0] == "custom":
            custom_lineage = st.text_input("Enter custom BUSCO lineage (e.g. fungi_odb10):", value="")
            lineage_code = custom_lineage.strip()
        else:
            lineage_code = busco_lineage[0]
        # BUSCO CPU selection (placed inside include_busco)
        max_cpus = os.cpu_count() or 4
        busco_cpus = st.slider("BUSCO CPUs", min_value=1, max_value=max_cpus, value=min(4, max_cpus))

    # Save uploaded files locally into session-specific temp directory
    # Ensure session_temp exists (fallback to environment var or 'temp_sessions')
    session_temp = st.session_state.get('session_temp') or os.environ.get('SESSION_TEMP_DIR')
    
    os.makedirs(os.path.dirname(session_temp), exist_ok=True)
    os.makedirs(session_temp, exist_ok=True)
    file_paths = file_handlers.save_uploaded_files(uploaded_files, directory_path=session_temp)
    st.session_state.temp_files_created = True
    st.session_state.session_temp = session_temp

    # Process uploaded files
    directory = st.session_state.get('session_temp') or os.environ.get('SESSION_TEMP_DIR')

    # Import here to avoid loading at startup
    from genome_analyzer import process_directory

    # If BUSCO analysis requested, require the user to confirm lineage and click Run
    df = None
    if include_busco:
        st.info("Choose a lineage and click 'Run analysis' to start BUSCO and other analyses.")
        run_now = st.button("Run analysis")
        if run_now:
            # Run per-file so we can show progress
            from genome_analyzer import compute_stats_with_busco
            with st.spinner("Analyzing genomes (including BUSCO)..."):
                results = []
                total = len(file_paths)
                progress = st.progress(0)
                status = st.empty()
                for idx, fp in enumerate(file_paths):
                    fname = os.path.basename(fp)
                    status.text(f"Running BUSCO for {fname} ({idx+1}/{total})")
                    stats = compute_stats_with_busco(fp, lineage=lineage_code, cpus=busco_cpus)
                    # ensure File column is present
                    if stats is None:
                        stats = {"File": fname.rsplit('.', 1)[0]}
                    else:
                        stats["File"] = fname.rsplit('.', 1)[0]
                    results.append(stats)
                    progress.progress(int(((idx+1) / total) * 100))
                # build dataframe from results and cache in session
                df = pd.DataFrame(results)
                st.session_state.qc_df = df
                # Mark that files have been processed for session cleanup
                st.session_state.files_processed = True
    else:
        # Non-BUSCO path can run immediately
        with st.spinner("Analyzing genomes..."):
            df = process_directory(directory, include_busco=False)
            # cache non-BUSCO results as well
            st.session_state.qc_df = df
            st.session_state.files_processed = True
    
    # Display results only when available
    # If we don't have a freshly computed df, check cached session results
    if df is None and st.session_state.get('qc_df') is not None:
        df = st.session_state.get('qc_df')

    if df is None:
        st.info("No analysis results yet. If you selected BUSCO, choose a lineage and click 'Run analysis'.")
    else:
        # Select columns to display
        if include_busco and any(col.startswith('BUSCO_') for col in df.columns):
            # Include BUSCO columns if available
            display_columns = ["File", "Total length (bp)", "Num contigs", "N50", "L90", "GC%", "BUSCO_Complete", "BUSCO_Single", "BUSCO_Fragmented", "BUSCO_Missing", "BUSCO_Total"]
            # Filter to only include columns that exist in the dataframe
            display_columns = [col for col in display_columns if col in df.columns]
            df_display = df[display_columns]
        else:
            # Basic columns only
            df_display = df[["File", "Total length (bp)", "Num contigs", "N50", "L90", "GC%"]]

        st.dataframe(df_display, hide_index=True)

        # Add explanations for metrics in expandable section
        with st.expander("What do these metrics mean?"):
            st.write("- **N50**: The length of the shortest contig at 50% of the total genome length when contigs are ordered by size")
            st.write("- **L90**: The number of contigs needed to cover 90% of the total genome length")
            if include_busco and any(col.startswith('BUSCO_') for col in df.columns):
                st.write("---")
                st.write("**BUSCO Metrics (Genome Completeness):**")
                st.write("- **BUSCO_Complete**: Percentage of complete BUSCO genes found (single + duplicated)")
                st.write("- **BUSCO_Fragmented**: Percentage of fragmented BUSCO genes")
                st.write("- **BUSCO_Missing**: Percentage of missing BUSCO genes")
                st.write("- **BUSCO_Total**: The total number of BUSCO genes searched for in the selected lineage")
                st.write("*Note: Higher completeness percentages indicate better genome quality*")

        # Export functionality
        if not df.empty:
            #st.write("---")
            st.write("")
            st.write("**Export Statistics Table**")
            
            # Custom filename input
            default_filename = "GenomeCheck_results"
            custom_filename = st.text_input(
                "Export filename (without extension):",
                value=default_filename,
                help="Enter a custom filename for your export. The appropriate extension (.csv or .xlsx) will be added automatically."
            )
            
            # Ensure filename is not empty
            if not custom_filename.strip():
                custom_filename = default_filename
            
            # Create two columns for the download buttons
            col1, col2 = st.columns(2)
            
            # CSV export
            csv_data = df_display.to_csv(index=False)
            with col1:
                st.download_button(
                    label="Download as CSV",
                    data=csv_data,
                    file_name=f"{custom_filename}.csv",
                    mime="text/csv"
                )
            
            # Excel export
            buffer = io.BytesIO()
            with pd.ExcelWriter(buffer, engine='openpyxl') as writer:
                df_display.to_excel(writer, index=False, sheet_name='GenomeCheck Results')
            excel_data = buffer.getvalue()
            
            with col2:
                st.download_button(
                    label="Download as Excel",
                    data=excel_data,
                    file_name=f"{custom_filename}.xlsx",
                    mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
                )

    #fastANI section
    st.write("---")
    st.subheader("All-vs-All fastANI Comparison & Similarity Visualization")

    if uploaded_files and len(file_paths) >= 2:
        # Create button for running analysis
        run_button = st.button("Run all-vs-all fastANI analysis")
        
        # Try to restore cached results from disk if they're not in session state
        if not st.session_state.get('fastani_done'):
            outputs_dir = os.path.join(st.session_state.get('session_temp') or os.environ.get('SESSION_TEMP_DIR', './temp/'), 'outputs')
            if os.path.exists(outputs_dir):
                # Restore heatmap bytes
                heatmap_bytes = {}
                tree_bytes = {}
                for fmt in ["png", "pdf", "svg", "jpeg"]:
                    # Restore heatmap
                    heatmap_path = os.path.join(outputs_dir, f"ani_heatmap.{fmt}")
                    if os.path.exists(heatmap_path):
                        with open(heatmap_path, 'rb') as f:
                            heatmap_bytes[fmt] = f.read()
                    
                    # Restore tree
                    tree_path = os.path.join(outputs_dir, f"nj_tree.{fmt}")
                    if os.path.exists(tree_path):
                        with open(tree_path, 'rb') as f:
                            tree_bytes[fmt] = f.read()
                
                # Restore tree newick text if it exists
                newick_path = os.path.join(outputs_dir, 'tree_newick.txt')
                if os.path.exists(newick_path):
                    with open(newick_path, 'r') as f:
                        st.session_state['tree_newick'] = f.read()

                # Restore genome names
                st.session_state.genome_names = [os.path.basename(fp) for fp in file_paths]

                if heatmap_bytes:
                    st.session_state['heatmap_bytes'] = heatmap_bytes
                    st.session_state['heatmap_preview'] = heatmap_bytes.get('png')
                    st.session_state['tree_bytes'] = tree_bytes
                    st.session_state['fastani_done'] = True

        # If we have cached results and the user did NOT press Run now, render cached outputs and downloads
        if st.session_state.get('fastani_done') and not run_button:
            st.write("**ANI Similarity Heatmap and Dendrogram (cached)**")
            # Show cached preview if available
            if st.session_state.get('heatmap_preview'):
                st.image(st.session_state.heatmap_preview)
            else:
                st.write("(No heatmap preview available)")

            # Show cached download buttons in one row
            buf_formats = [
                ("png", "image/png"),
                ("pdf", "application/pdf"),
                ("svg", "image/svg+xml"),
                ("jpeg", "image/jpeg"),
            ]
            cols = st.columns(len(buf_formats))
            for i, (fmt, mime) in enumerate(buf_formats):
                with cols[i]:
                    data = st.session_state.get('heatmap_bytes', {}).get(fmt) if st.session_state.get('heatmap_bytes') else None
                    if data:
                        st.download_button(label=f"{fmt.upper()}", data=data, file_name=f"ani_heatmap.{fmt}", mime=mime, key=f"heatmap_download_{fmt}")
                    else:
                        st.write("-")

            # Tree section from cache
            st.write("### Neighbor-Joining Tree")
            if len(st.session_state.get('genome_names', [])) < 3:
                st.info("Tree not generated: need at least 3 genomes to build a neighbor-joining tree.")
            elif st.session_state.get('tree_bytes'):
                tree_preview = st.session_state['tree_bytes'].get('png')
                if tree_preview:
                    st.image(tree_preview)
                else:
                    st.write("(No tree preview available)")

                tree_cols = st.columns(len(buf_formats))
                for i, (fmt, mime) in enumerate(buf_formats):
                    with tree_cols[i]:
                        data = st.session_state['tree_bytes'].get(fmt)
                        if data:
                            st.download_button(label=f"{fmt.upper()}", data=data, file_name=f"nj_tree.{fmt}", mime=mime, key=f"tree_download_{fmt}")
                        else:
                            st.write("-")

                if len(st.session_state.get('genome_names', [])) >= 3:
                    with st.expander("Newick format text (cached)"):
                        st.code(st.session_state.get('tree_newick', ''))

        if run_button:
            # Clear previous cached data before new analysis
            st.session_state.pop('heatmap_bytes', None)
            st.session_state.pop('heatmap_preview', None)
            st.session_state.pop('tree_bytes', None)
            st.session_state.pop('tree_newick', None)
            st.session_state['fastani_done'] = False
            
            with st.spinner("Running fastANI comparisons..."):
                ani_matrix, genome_names, fastani_errors = all_vs_all_fastani(file_paths)
            if fastani_errors:
                st.warning("Some pairwise fastANI comparisons failed. See details below.")
                for (i, j), msg in fastani_errors.items():
                    st.write(f"Comparison failed: {os.path.basename(file_paths[i])} vs {os.path.basename(file_paths[j])}: {msg}")
            st.success("fastANI analysis complete!")
            # Replace NaNs with 0 for clustering
            ani_matrix = np.nan_to_num(ani_matrix, nan=0.0)
            # Show clustermap (dendrogram + heatmap)
            st.write("**ANI Similarity Heatmap and Dendrogram**")
            cg = sns.clustermap(
                ani_matrix,
                metric="euclidean",
                method="average",
                cmap="viridis",
                annot=True,
                fmt=".2f",
                xticklabels=genome_names,
                yticklabels=genome_names,
                figsize=(10, 8)
            )
            
            # Save figure as bytes for preview
            preview_buf = io.BytesIO()
            cg.figure.savefig(preview_buf, format='png', bbox_inches='tight')
            preview_buf.seek(0)
            st.session_state.heatmap_preview = preview_buf.getvalue()
            
            # Display the heatmap
            st.image(st.session_state.heatmap_preview)

            # Export DPI selector for raster formats
            dpi = st.selectbox("Export resolution (DPI) for PNG/JPEG", options=[72, 150, 300], index=1)

            # Prepare/download formats and cache bytes in session_state to avoid recomputation on widget reruns
            buf_formats = [
                ("png", "image/png"),
                ("pdf", "application/pdf"),
                ("svg", "image/svg+xml"),
                ("jpeg", "image/jpeg"),
            ]

            # Save all formats and update session state
            heatmap_bytes = {}
            outputs_dir = os.path.join(st.session_state.get('session_temp') or os.environ.get('SESSION_TEMP_DIR', './temp/'), 'outputs')
            os.makedirs(outputs_dir, exist_ok=True)
            
            for fmt, mime in buf_formats:
                try:
                    buf = io.BytesIO()
                    if fmt in ("png", "jpeg"):
                        cg.figure.savefig(buf, format=fmt, bbox_inches='tight', dpi=dpi)
                    else:
                        cg.figure.savefig(buf, format=fmt, bbox_inches='tight')
                    buf.seek(0)
                    bytes_data = buf.getvalue()
                    heatmap_bytes[fmt] = bytes_data
                    
                    # Save to disk for persistence
                    with open(os.path.join(outputs_dir, f"ani_heatmap.{fmt}"), 'wb') as f:
                        f.write(bytes_data)
                except Exception as e:
                    heatmap_bytes[fmt] = None
                    st.warning(f"Could not create {fmt} for heatmap: {e}")

            # Update session state
            st.session_state['heatmap_bytes'] = heatmap_bytes
            st.session_state['fastani_done'] = True

            # Show downloads in one row using cached bytes
            cols = st.columns(len(buf_formats))
            for i, (fmt, mime) in enumerate(buf_formats):
                with cols[i]:
                    data = st.session_state['heatmap_bytes'].get(fmt)
                    if data:
                        st.download_button(label=f"{fmt.upper()}", data=data, file_name=f"ani_heatmap.{fmt}", mime=mime, key=f"cached_heatmap_download_{fmt}")
                    else:
                        st.write("-")
            
            # Show neighbor-joining tree (rendered with Biopython/Matplotlib)
            st.write("### Neighbor-Joining Tree")
            if len(genome_names) < 3:
                st.info("Tree not generated: need at least 3 genomes to build a neighbor-joining tree.")
            else:
                distance_matrix = 1 - (ani_matrix / 100.0)
                tree_ascii, tree_newick = neighbor_joining_tree(distance_matrix, genome_names)
                if tree_newick:
                    try:
                        tree_obj = Phylo.read(StringIO(tree_newick), "newick")
                        tree_obj.ladderize()
                        # size depends on number of taxa
                        height = max(4, len(genome_names) * 0.3)
                        fig = plt.figure(figsize=(6, height))
                        ax = fig.add_subplot(1, 1, 1)
                        Phylo.draw(tree_obj, do_show=False, axes=ax, show_confidence=False, label_func=lambda n: n.name)
                        # remove axes, ticks and spines for a cleaner tree-only view
                        ax.set_axis_off()
                        for spine in getattr(ax, 'spines', {}).values():
                            spine.set_visible(False)
                        plt.tight_layout()
                        # Show cached tree preview if available
                        tree_preview = None
                        if st.session_state.get('tree_bytes'):
                            tree_preview = st.session_state['tree_bytes'].get('png')
                        if tree_preview:
                            st.image(tree_preview)
                        else:
                            st.pyplot(fig)

                        # Cache tree images and show download buttons from cache
                        tree_bytes = {}
                        for fmt, mime in buf_formats:
                            try:
                                buf = io.BytesIO()
                                if fmt in ("png", "jpeg"):
                                    fig.savefig(buf, format=fmt, bbox_inches='tight', dpi=dpi)
                                else:
                                    fig.savefig(buf, format=fmt, bbox_inches='tight')
                                buf.seek(0)
                                tree_bytes[fmt] = buf.getvalue()
                            except Exception as e:
                                tree_bytes[fmt] = None
                                st.warning(f"Could not create {fmt} for tree: {e}")
                        # Save tree newick to disk for persistence
                        with open(os.path.join(outputs_dir, 'tree_newick.txt'), 'w') as f:
                            f.write(tree_newick)
                        for fmt in buf_formats:
                            # Save to disk
                            with open(os.path.join(outputs_dir, f"nj_tree.{fmt[0]}"), 'wb') as f:
                                f.write(tree_bytes[fmt[0]])
                        st.session_state['tree_bytes'] = tree_bytes 
                        tree_cols = st.columns(len(buf_formats))
                        for i, (fmt, mime) in enumerate(buf_formats):
                            with tree_cols[i]:
                                data = st.session_state['tree_bytes'].get(fmt)
                                if data:
                                    st.download_button(label=f"{fmt.upper()}", data=data, file_name=f"nj_tree.{fmt}", mime=mime, key=f"cached_tree_download_{fmt}")
                                else:
                                    st.write("-")

                        # Mark that fastANI results are available in session state
                        st.session_state['fastani_done'] = True
                        st.session_state['ani_matrix'] = ani_matrix
                        st.session_state['genome_names'] = genome_names
                    except Exception as e:
                        st.write("Failed to render tree with Biopython; falling back to ASCII output")
                        st.code(tree_ascii)
                else:
                    st.write("Neighbor-joining tree not available.")

            if len(genome_names) >= 3:
                with st.expander("Newick format text"):
                    st.code(tree_newick or "")

    # Session-based cleanup
    if 'files_processed' not in st.session_state:
        st.session_state.files_processed = False
        st.session_state.temp_files_created = False
    
    # Clean up on page refresh if files were created in this session
    if st.session_state.temp_files_created and not uploaded_files:
        # If no files are uploaded but we had temp files, clean up this session's temp dir
        file_handlers.cleanup_temp_directory()
        st.session_state.temp_files_created = False
        st.session_state.files_processed = False
    
    # Manual cleanup option
    if st.button("Clean up temporary files"):
        file_handlers.cleanup_temp_directory()
        st.success("Temporary files cleaned up!")
        st.rerun()
    
    # Register automatic cleanup (runs when Python process exits)
    file_handlers.register_cleanup_on_exit(session_temp)
