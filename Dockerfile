# Start from the minimal conda-forge image
FROM condaforge/mambaforge:latest

# Set working directory
WORKDIR /app

# Copy your repo contents
COPY . .

# Create the environment
#Include BUSCO or not
ARG INCLUDE_BUSCO=false
RUN if [ "$INCLUDE_BUSCO" = "true" ]; then \
      echo "Installing BUSCO and dependencies..."; \
      mamba install -y -c conda-forge -c bioconda busco=6.0.0; \
    else \
      echo "Skipping BUSCO to keep image small"; \
    fi

RUN pip install --no-cache-dir streamlit biopython openpyxl scikit-bio seaborn
RUN mamba install -c conda-forge -c bioconda fastani

# Expose Streamlit port
EXPOSE 8501

# Default command to run your Streamlit app
CMD ["streamlit", "run", "streamlit_app.py", "--server.port", "8501", "--server.address", "0.0.0.0"]
