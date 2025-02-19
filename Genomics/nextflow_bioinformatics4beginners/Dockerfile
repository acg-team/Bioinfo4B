# Use Miniconda3 base image
FROM continuumio/miniconda3:latest

# Set environment variables
ENV PATH="/opt/conda/bin:$PATH" \
    DEBIAN_FRONTEND=noninteractive

# Install Mamba for faster package management
RUN conda install -y -c conda-forge mamba && \
    mamba update -y conda && \
    mamba install -y -c conda-forge -c bioconda -c defaults \
        python=3.8 \
        gangstr \
        bwa-mem2 \
        cyvcf2 \
        figtree \
        samtools \
        phylip \
        numpy \
        pandas \
        seaborn \
        biopython \
        notebook \
        fastx_toolkit \
        fastqc \
        python-wget \
        bowtie2 \
        bcftools \
        clustalo \
        cnvkit && \
    conda clean -afy

# Set default shell to bash and activate conda automatically
SHELL ["/bin/bash", "-c"]

# Ensure conda environment is always activated
RUN echo "source activate base" >> ~/.bashrc
ENV PATH="/opt/conda/bin:$PATH"

# Define default command
CMD ["bash"]
