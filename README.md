# Genome QC App

An easy-to-use Streamlit application for basic genome quality control (QC) and similarity analysis. It computes common assembly metrics (N50, L90, GC%), optionally runs BUSCO for genome completeness, and can compute pairwise ANI with fastANI and show a clustered heatmap + tree.

This README is written for users who may be new to the command line and bioinformatics tools.

## Quick overview

- Upload one or more genome FASTA files (.fasta, .fa, .fna)
- Optionally run BUSCO (select lineage and click "Run analysis")
- Run all-vs-all fastANI to produce a clustered ANI heatmap and neighbor-joining tree
- Export results as CSV or Excel

## Contents of this repository

- `streamlit_app.py` — Streamlit user interface
- `genome_analyzer.py` — Core analysis functions (statistics, BUSCO wrapper, fastANI helper)
- `utils/file_handlers.py` — temporary file handling and cleanup helpers
- `.streamlit/config.toml` — Streamlit configuration used by the app
- `requirements.txt` — Python dependencies (see notes below)
- `DEPLOYMENT.md` — deployment instructions for Docker / Streamlit Cloud

## Installation (local development)

These steps assume you have Python 3.8+ installed. If you're unfamiliar with Python environments, the easiest path on macOS is to use Miniforge/Miniconda.

1. (Optional but recommended) Create and activate a Python environment (conda example):

```bash
conda create -n genome_qc python=3.10 -y
conda activate genome_qc
```

2. Install Python requirements:

```bash
pip install -r requirements.txt
```

3. Install external tools used by the app (not Python packages):

- BUSCO (optional, for completeness checks). Recommended install with conda:

```bash
conda install -c bioconda busco
```

- fastANI (optional, for ANI comparisons):

```bash
conda install -c bioconda fastani
```

4. Run the Streamlit app locally:

```bash
streamlit run streamlit_app.py
```

Open the URL printed by Streamlit (usually http://localhost:8501) in your browser.

If you prefer Docker or Streamlit Cloud, see `DEPLOYMENT.md`.

## How to use the app (step-by-step)

1. Open the app in your browser.
2. Use the "Upload genome FASTA files" control to upload one or more assemblies.
3. (Optional) Tick "Include BUSCO analysis" if you want completeness estimates.
   - If you do, choose a lineage from the dropdown or select "Custom" and paste your lineage name.
   - Choose the number of CPUs with the slider (default: 4). This controls BUSCO parallelism and speeds up runs on multi-core machines.
   - Click the "Run analysis" button to start BUSCO; a progress bar and per-file status updates will show progress.
4. After analysis, view the results table and explanations.
5. (Optional) Run the "Run all-vs-all fastANI analysis" button to compute pairwise ANI across uploaded genomes and generate a clustered heatmap + tree.
6. Export results as CSV or Excel using the download buttons.

## Notes about BUSCO

- BUSCO must be installed separately (conda bioconda channel recommended).
- BUSCO can be slow for large genomes; increase CPUs to speed it up, but be mindful of memory.
- BUSCO expects lineage datasets to be available. The app will use BUSCO's configuration/location; you can pre-download lineages to avoid network downloads.
- If you provide a "Custom" lineage, make sure that lineage name is valid and available to BUSCO on the machine where the app runs.

## Notes about fastANI

- fastANI is executed as an external program (subprocess) and does not import numpy/pandas internals — this avoids package version conflicts.
- Ensure `fastANI` is on PATH or installed in the same environment used to run Streamlit.

## Troubleshooting

- "BUSCO not found" errors: install BUSCO in the same environment you use to run Streamlit.
- "Data must be symmetric and cannot contain NaNs" (tree building): some fastANI comparisons may fail; re-run or inspect logs. The app replaces NaNs with zero for clustering.
- If the app cannot find `fastANI` or `busco`, test from your terminal:

```bash
busco --version
fastANI --version
```

If these fail, install the tool or adjust your PATH.

## Advanced tips

- For big batches of genomes, consider running BUSCO on a cluster or machine with more RAM/CPUs and importing results into the app.
- To speed up BUSCO, install DIAMOND (bioconda) so BUSCO will use it for sequence searches when appropriate.
- To create publication-quality trees, you can export the Newick string from the app and use dedicated tree viewers (iTOL, ETE3).

## Example workflow (for novices)

1. Install Miniforge/conda and create the `genome_qc` environment (see Installation).
2. Install prerequisites and tools.
3. Launch the app: `streamlit run streamlit_app.py`.
4. Upload 2–6 small bacterial assemblies (example FASTA files can be downloaded from NCBI).
5. Tick "Include BUSCO", choose "Bacteria (odb12)", set CPUs to 4, click "Run analysis".
6. When BUSCO finishes, click "Run all-vs-all fastANI analysis" and inspect the heatmap + tree.

## Contributing and development

Contributions welcome — please open issues or pull requests. If you change behavior, add/update tests where appropriate.

## License

This project is released under the MIT License - see `LICENSE` for details.
