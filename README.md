# Genome QC App

An easy-to-use Streamlit application for basic genome quality control (QC) and similarity analysis. The app will compute common assembly metrics (N50, L90, GC%), optionally runs [BUSCO](https://busco.ezlab.org/) for genome completeness, and can estimate pairwise Average Nucleotide Identity (ANI) with [fastANI](https://www.nature.com/articles/s41467-018-07641-9) before showing a clustered heatmap and tree.

## Quick overview

- Upload one or more genome FASTA files (.fasta, .fa, .fna)
- Optionally run BUSCO (select or provide lineage to use as the source of BUSCOs and click "Run analysis")
- Optionally run all-vs-all fastANI to produce a clustered ANI heatmap and neighbor-joining tree
- Export tabular results as a CSV or Excel file and image results as PNG, PDF, SVG, or JPEG

## Streamlit Community Cloud

Use the hosted app here: [genomeqc.streamlit.app](https://genomeqc.streamlit.app/).

Notes:
- BUSCO is disabled on the hosted Cloud app. For BUSCO features, run locally or use the BUSCO-enabled Docker image (see below).

## Run with Docker

*Running the app via Docker or a local installation is recommended if BUSCO analysis or more compute resources are needed (e.g. genome upload for large FASTAs or fastANI analyses are running slow on Streamlit Community Cloud)*

### Quick start (recommended)

**Basic version (no BUSCO):**
```bash
docker run --rm -p 8501:8501 meadm/genome-qc:latest
```

**BUSCO-enabled version:**
```bash
docker run --rm -p 8501:8501 -v "$PWD:/app/work" meadm/genome-qc:busco
```

Then open the URL the command line displays (it should be something like http://0.0.0.0:8501) in your browser.

## Local Installation and Use

These steps assume you have conda/mamba available (e.g., via [Miniforge3](https://github.com/conda-forge/miniforge)).

0. Clone this repository and enter the project directory:

```bash
git clone https://github.com/meadm/genome_QC.git
cd genome_QC
```

1. Create the environment from the `environment.yml` file and activate it:

```bash
mamba env create -f environment.yml  # or: conda env create -f environment.yml
conda activate genome_qc
```

2. Install Streamlit (if not already present in the environment):

```bash
pip install streamlit
```

3. BUSCO is NOT included in the distributed `environment.yml` to keep installs lighter. If you need BUSCO locally, install it separately:

```bash
conda install -c conda-forge -c bioconda busco=6.0.0
```

4. Run the app:

```bash
streamlit run streamlit_app.py
```

Open the URL printed by Streamlit (usually http://localhost:8501) in your web browser.

## How to use the app (step-by-step)

1. Open the app in your browser.
2. Use the "Upload genome FASTA files" section to upload one or more assemblies.
3. (Optional and disabled on Streamlit Community Cloud) Tick "Include BUSCO analysis" if you want completeness estimates.
   - If you do, choose a lineage from the dropdown or select "Custom" and paste your lineage name.
   - Choose the number of CPUs with the slider (default: 4). This controls BUSCO parallelism and speeds up runs on multi-core machines.
   - Click the "Run analysis" button to start BUSCO; a progress bar and per-file status updates will show progress.
   - Export tabular results as CSV or Excel using the download buttons.
4. (Optional) Run the "Run all-vs-all fastANI analysis" button to compute pairwise ANI across uploaded genomes and generate a clustered heatmap + tree.
   - Note that fastANI is *estimating* whole genome simlilarity so the output trees may not be phylogenetically correct.
   - Download result images as PNG, PDF, SVG, or JPEG

## Notes about BUSCO

- BUSCO must be installed separately (conda/bioconda recommended) when running locally.
- BUSCO can be slow for large (> 15Mb) genomes; increase CPUs to speed it up, but be mindful of local memory.
- BUSCO expects lineage datasets to be available. The app will use BUSCO's configuration/location; you can pre-download lineages to avoid network downloads.
- If you provide a "Custom" lineage, make sure that lineage name is valid and available to BUSCO on the machine where the app runs.

## For developers

To build Docker images locally:
(Assuming you have cloned the repo and are in the project directory)

```bash
# Lean image (default; no BUSCO)
docker build -t genome-qc:latest .

# BUSCO-enabled image
docker build --build-arg INCLUDE_BUSCO=true -t genome-qc:busco .
```
Then run using the `docker run` commands above, but omit the `meadm` prefix to the image name.

## Contributing and development

Contributions, feedback, and feature requests are welcome! Please:

- For developers, create a feature branch from `main` (e.g., `feature/better-export`, `fix/tree-guard`).
- Add clear descriptions of the changes and rationales to Pull Requests.
- Prefer "Squash and merge" to keep history tidy.
- For larger changes or requests, open an issue first to discuss scope and approach.

## License

This project is released under the MIT License - see `LICENSE` for details.
