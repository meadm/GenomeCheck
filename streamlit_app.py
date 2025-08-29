import streamlit as st
import pandas as pd
from utils import file_handlers

st.title("Genome QC and Similarity App")

# -------------------------
# 1. File Upload
# -------------------------
uploaded_files = st.file_uploader(
    "Upload genome FASTA files", type=["fasta", "fa", "fna"], accept_multiple_files=True
)

if uploaded_files:
    st.write(f"{len(uploaded_files)} genomes uploaded")

    # BUSCO options
    st.write("---")
    st.subheader("Analysis Options")
    include_busco = st.checkbox(
        "Include BUSCO analysis (genome completeness assessment)",
        value=False,
        help="BUSCO assesses genome completeness by searching for universal single-copy orthologs. This may take several minutes per genome."
    )
    
    if include_busco:
        st.info("⚠️ **Note**: BUSCO analysis can take 5-30 minutes per genome depending on size and complexity.")
        
        # Lineage selection
        busco_lineage = st.selectbox(
            "Select organism lineage for BUSCO analysis:",
            options=[
                ("bacteria_odb12", "Bacteria (odb12)"),
                ("eukaryota_odb12", "Eukaryota (odb12)"),
                ("archaea_odb12", "Archaea (odb12)")
            ],
            format_func=lambda x: x[1],
            help="Choose the appropriate lineage for your organisms. Bacteria for prokaryotes, Eukaryota for eukaryotes, Archaea for archaea."
        )
        
        # Extract the lineage code
        lineage_code = busco_lineage[0]

    # Save uploaded files locally
    file_handlers.create_temp_directory()
    file_paths = file_handlers.save_uploaded_files(uploaded_files)
    st.session_state.temp_files_created = True

    # Process uploaded files
    directory = './temp/'
    
    with st.spinner("Analyzing genomes..."):
        # Import here to avoid loading at startup
        from genome_analyzer import process_directory
        
        if include_busco:
            df = process_directory(directory, include_busco=True, busco_lineage=lineage_code)
        else:
            df = process_directory(directory, include_busco=False)
        
        # Mark that files have been processed for session cleanup
        st.session_state.files_processed = True
    
    # Select columns to display
    if include_busco and any(col.startswith('BUSCO_') for col in df.columns):
        # Include BUSCO columns if available
        display_columns = ["File", "Total length", "Num contigs", "N50", "L90", "GC%", "BUSCO_Complete", "BUSCO_Single", "BUSCO_Fragmented", "BUSCO_Missing", "BUSCO_Total"]
        # Filter to only include columns that exist in the dataframe
        display_columns = [col for col in display_columns if col in df.columns]
        df = df[display_columns]
    else:
        # Basic columns only
        df = df[["File", "Total length", "Num contigs", "N50", "L90", "GC%"]]

    st.dataframe(df, hide_index=True)
    
    # Add explanations for metrics in expandable section
    with st.expander("What do these metrics mean?"):
        st.write("- **N50**: The length of the shortest contig at 50% of the total genome length when contigs are ordered by size")
        st.write("- **L90**: The number of contigs needed to cover 90% of the total genome length")
        if include_busco and any(col.startswith('BUSCO_') for col in df.columns):
            st.write("---")
            st.write("**BUSCO Metrics (Genome Completeness):**")
            st.write("- **BUSCO_Complete**: Percentage of complete BUSCO genes found (single + duplicated)")
            st.write("- **BUSCO_Single**: Percentage of complete single-copy BUSCO genes")
            st.write("- **BUSCO_Fragmented**: Percentage of fragmented BUSCO genes")
            st.write("- **BUSCO_Missing**: Percentage of missing BUSCO genes")
            st.write("- **BUSCO_Total**: The total number of BUSCO genes searched for in the selected lineage")
            st.write("*Note: Higher completeness percentages indicate better genome quality*")

    # Export functionality
    if not df.empty:
        st.write("---")
        st.subheader("Export Results")
        
        # Custom filename input
        default_filename = "genome_qc_results"
        custom_filename = st.text_input(
            "Export filename (without extension):",
            value=default_filename,
            help="Enter a custom filename for your export. The appropriate extension (.csv or .xlsx) will be added automatically."
        )
        
        # Ensure filename is not empty
        if not custom_filename.strip():
            custom_filename = default_filename
        
        # CSV export
        csv_data = df.to_csv(index=False)
        st.download_button(
            label="Download as CSV",
            data=csv_data,
            file_name=f"{custom_filename}.csv",
            mime="text/csv"
        )
        
        # Excel export
        import io
        buffer = io.BytesIO()
        with pd.ExcelWriter(buffer, engine='openpyxl') as writer:
            df.to_excel(writer, index=False, sheet_name='Genome QC Results')
        excel_data = buffer.getvalue()
        
        st.download_button(
            label="Download as Excel",
            data=excel_data,
            file_name=f"{custom_filename}.xlsx",
            mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
        )

    # Session-based cleanup
    if 'files_processed' not in st.session_state:
        st.session_state.files_processed = False
        st.session_state.temp_files_created = False
    
    # Clean up on page refresh if files were created in this session
    if st.session_state.temp_files_created and not uploaded_files:
        # If no files are uploaded but we had temp files, clean up
        file_handlers.cleanup_temp_directory()
        st.session_state.temp_files_created = False
        st.session_state.files_processed = False
    
    # Manual cleanup option
    if st.button("Clean up temporary files"):
        file_handlers.cleanup_temp_directory()
        st.success("Temporary files cleaned up!")
        st.rerun()
    
    # Register automatic cleanup (runs when Python process exits)
    file_handlers.register_cleanup_on_exit()