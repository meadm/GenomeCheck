import streamlit as st
import pandas as pd
from Bio import SeqIO
import os
import shutil
import atexit

st.title("Genome QC and Similarity App")

# -------------------------
# 1. File Upload
# -------------------------
uploaded_files = st.file_uploader(
    "Upload genome FASTA files", type=["fasta", "fa", "fna"], accept_multiple_files=True
)

if uploaded_files:
    st.write(f"{len(uploaded_files)} genomes uploaded")

    # Save uploaded files locally
    try:
        os.makedirs('./temp/', exist_ok=True)
    except OSError as e:
        print(f"Error creating directory: {e}")

    file_paths = []
    for f in uploaded_files:
        path = f"./temp/temp_{f.name}"
        with open(path, "wb") as out:
            out.write(f.read())
        file_paths.append(path)

    def compute_stats(fasta_file):
        lengths = []
        gc_count = 0
        total_bases = 0

        for record in SeqIO.parse(fasta_file, "fasta"):
            seq = record.seq.upper()
            seq_len = len(seq)
            lengths.append(seq_len)
            gc_count += seq.count("G") + seq.count("C")
            total_bases += seq_len

        if not lengths:
            return {"Total length": 0, "Num contigs": 0, "N50": 0, "L90": 0, "GC%": 0}

        lengths.sort(reverse=True)
        total_length = sum(lengths)
        num_contigs = len(lengths)
        gc_percent = (gc_count / total_bases) * 100 if total_bases > 0 else 0

        # N50
        half_length = total_length / 2
        cumulative = 0
        n50 = 0
        for length in lengths:
            cumulative += length
            if cumulative >= half_length:
                n50 = length
                break

        # L90
        l90_target = total_length * 0.9
        cumulative = 0
        l90 = 0
        for length in lengths:
            cumulative += length
            l90 += 1
            if cumulative >= l90_target:
                break

        return {
            "Total length": total_length,
            "Num contigs": num_contigs,
            "N50": n50,
            "L90": l90,
            "GC%": round(gc_percent, 2)
        }

    def process_directory(directory):
        results = []
        for filename in os.listdir(directory):
            if filename.endswith((".fasta", ".fa", ".fna")):
                filepath = os.path.join(directory, filename)
                stats = compute_stats(filepath)
                stats["File"] = filename
                results.append(stats)
        return pd.DataFrame(results)

    # Example usage
    directory = './temp/'
    df = process_directory(directory)
    df = df[["File", "Total length", "Num contigs", "N50", "L90", "GC%"]]  # reorder columns

    st.dataframe(df)

    def cleanup():
        shutil.rmtree('./temp/', ignore_errors=True)

    atexit.register(cleanup)