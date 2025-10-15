# Start from the minimal conda-forge image
FROM condaforge/mambaforge:latest

# Set working directory
WORKDIR /app

# Copy your repo contents
COPY . .

# Create the environment
#RUN conda install -c conda-forge -c bioconda busco=6.0.0
RUN conda install -c conda-forge -c bioconda compleasm
RUN pip install streamlit
RUN pip install biopython
RUN pip install openpyxl
RUN conda install -c conda-forge -c bioconda fastani
RUN pip install scikit-bio seaborn

# Expose Streamlit port
EXPOSE 8501

# Default command to run your Streamlit app
CMD ["streamlit", "run", "streamlit_app.py", "--server.port", "8501", "--server.address", "0.0.0.0"]
