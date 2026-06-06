
# Do this if you suffer from low disk space
#   DOCKER_TMPDIR=/mnt/data/docker-tmp

# 1. Base Image: NVIDIA's official CUDA toolkit (Development version includes compilers)
FROM nvidia/cuda:12.4.1-devel-ubuntu22.04

ARG JUPYTER_HASH
ENV JUPYTER_PASSWORD_HASH=$JUPYTER_HASH

# 2. Set environment variables to prevent interactive prompts
ENV DEBIAN_FRONTEND=noninteractive
ENV CONDA_DIR=/opt/conda
ENV PATH=$CONDA_DIR/bin:$PATH

# 3. Install system dependencies
RUN apt-get update && apt-get install -y \
    zsh \
    curl \
    git \
    nvtop \
    wget \
    bzip2 \
    ca-certificates \
    && rm -rf /var/lib/apt/lists/*

# 4. Install Oh My Zsh (Unattended)
RUN sh -c "$(curl -fsSL https://raw.githubusercontent.com/ohmyzsh/ohmyzsh/master/tools/install.sh)" "" --unattended

# 5. Change default shell to Zsh
RUN chsh -s $(which zsh)

# 6. Install Miniconda
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /tmp/miniconda.sh && \
    /bin/bash /tmp/miniconda.sh -b -p $CONDA_DIR && \
    rm /tmp/miniconda.sh && \
    $CONDA_DIR/bin/conda clean -afy

# 7. Initialize conda for zsh and bash
RUN conda init zsh && conda init bash

# 8. Create the conda environment with Python 3.13, RAPIDS 26.02, and JupyterLab
RUN conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main && \
    conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r && \
    conda create -y -n rapids-env -c rapidsai -c conda-forge -c nvidia \
    rapids=26.02 python=3.13 cuda-version=12.4 jupyterlab && \
    conda clean -afy

RUN /opt/conda/envs/rapids-env/bin/python -m ipykernel install --sys-prefix --name rapids-env --display-name "Python 3.13 (RAPIDS)"

# 9. Ensure the conda environment activates automatically when you SSH/exec into Zsh
RUN echo "conda activate rapids-env" >> ~/.zshrc

RUN apt-get update && apt-get install -y htop lsof && \
    rm -rf /var/lib/apt/lists/* && \
    conda install -y -n rapids-env -c conda-forge nvitop

RUN git clone https://github.com/perpk/methylation-analysis.git

# 10. Expose the Jupyter port
EXPOSE 8080

# 11. Launch Jupyter Lab within the RAPIDS environment when the container starts
CMD ["sh", "-c", "conda run --no-capture-output -n rapids-env jupyter lab --ip=0.0.0.0 --port=8080 --no-browser --allow-root --ServerApp.password=\"$JUPYTER_PASSWORD_HASH\""]