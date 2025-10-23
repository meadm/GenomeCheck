# Multi-stage build for smaller final images
# Stage 1: Build and install everything
FROM condaforge/mambaforge:latest AS builder

WORKDIR /app
COPY . .

# Install packages based on BUSCO flag
ARG INCLUDE_BUSCO=false
RUN if [ "$INCLUDE_BUSCO" = "true" ]; then \
      echo "Installing BUSCO and dependencies..."; \
      mamba install -y -c conda-forge -c bioconda busco=6.0.0 fastani scikit-bio seaborn; \
    else \
      echo "Installing lean dependencies..."; \
      mamba install -y -c conda-forge -c bioconda fastani scikit-bio seaborn; \
    fi

# Install Python packages
RUN pip install --no-cache-dir streamlit biopython openpyxl

# Clean up conda and pip caches
RUN mamba clean --all -y && pip cache purge

# Stage 2: Create lean runtime image
FROM condaforge/mambaforge:latest

WORKDIR /app

# Copy installed packages from builder stage
COPY --from=builder /opt/conda /opt/conda

# Copy application code
COPY . .

# Create non-root user for security
RUN useradd -m -u 1000 appuser && chown -R appuser:appuser /app
USER appuser

# Set environment variables
ENV BUSCO_DATA_DIR=/app/work/busco_downloads
ENV PATH=/opt/conda/bin:$PATH

# Expose Streamlit port
EXPOSE 8501

# Health check
HEALTHCHECK --interval=30s --timeout=5s --retries=3 \
  CMD wget -qO- http://localhost:8501/_stcore/health || exit 1

# Default command to run your Streamlit app
CMD ["streamlit", "run", "streamlit_app.py", "--server.port", "8501", "--server.address", "0.0.0.0"]
