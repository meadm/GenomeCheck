# Genome QC App

A comprehensive genome quality control and analysis application with support for both Docker and Streamlit Cloud deployment.

## Features

- **Basic Genome QC**: N50, L90, GC content, contig statistics
- **BUSCO Analysis**: Genome completeness assessment (optional)
- **Multiple Lineages**: Support for bacteria, eukaryota, and archaea
- **Export Options**: CSV and Excel export with custom filenames
- **User-Friendly**: Interactive interface with helpful explanations

## Quick Start

### Option 1: Streamlit Cloud (Recommended for testing)
Visit the deployed app: [Link to be added after deployment]

### Option 2: Docker (Recommended for production)
```bash
docker run -p 8501:8501 yourusername/genome-qc-app:latest
```

### Option 3: Local Development
```bash
pip install -r requirements.txt
streamlit run streamlit_app.py
```

## Deployment

See [DEPLOYMENT.md](DEPLOYMENT.md) for detailed deployment instructions for both Docker and Streamlit Cloud.

## Usage

1. Upload FASTA files (.fasta, .fa, .fna)
2. Choose analysis options (BUSCO analysis is optional)
3. Select appropriate lineage for BUSCO analysis
4. View results in the interactive table
5. Export results as CSV or Excel

## Architecture

- `streamlit_app.py`: Main UI and user interface
- `genome_analyzer.py`: Core analysis logic
- `utils/file_handlers.py`: File operations and cleanup
