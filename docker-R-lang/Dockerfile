# 1. Base Image: Rocker's RStudio image (comes with R, RStudio Server, and tidyverse)
FROM rocker/rstudio:latest

# 2. Configure RStudio Server for vast.ai
# Vast.ai's default web proxy looks for port 8080. We map RStudio to listen there.
RUN echo "www-port=8080" >> /etc/rstudio/rserver.conf

# Disable authentication so you can click the vast.ai "Connect" button and jump straight in.
# (Vast.ai's proxy link is already secure/authenticated to your account)
ENV DISABLE_AUTH=true

# 3. Install critical system dependencies
# When processing large multi-omics datasets or compiling heavy Bioconductor packages 
# from source, R requires specific underlying C++ libraries. Missing these is the #1 
# cause of R package installation failures.
RUN apt-get update && apt-get install -y \
    openssh-server \
    libcurl4-openssl-dev \
    libxml2-dev \
    libssl-dev \
    libfontconfig1-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libpng-dev \
    zlib1g-dev \
    cmake \
    htop \
    lsof \
    zsh \
    curl \
    git \
    nvtop \
    wget \
    bzip2 \
    ca-certificates \
    && rm -rf /var/lib/apt/lists/*

RUN sh -c "$(curl -fsSL https://raw.githubusercontent.com/ohmyzsh/ohmyzsh/master/tools/install.sh)" "" --unattended

RUN chsh -s $(which zsh)

# 4. Set the working directory
WORKDIR /workspace

# 5. Copy your custom R install script into the container
# Make sure "install_script.R" is in the same folder as this Dockerfile when building
COPY docker-R-lang/install_script.R /workspace/install_script.R

# 6. Execute the R script during the build phase
# This bakes all your downloaded libraries directly into the image layers
RUN Rscript /workspace/install_script.R

RUN git clone https://github.com/perpk/methylation-analysis.git

# Configure SSH to disable password logins and explicitly allow root access
RUN mkdir /var/run/sshd && \
    echo "PasswordAuthentication no" >> /etc/ssh/sshd_config && \
    echo "PermitRootLogin yes" >> /etc/ssh/sshd_config

# Tell Rocker's init system to handle the Vast.ai SSH key and start the daemon
RUN mkdir -p /etc/services.d/sshd && \
    printf '#!/bin/bash\n\
mkdir -p /root/.ssh\n\
touch /root/.ssh/authorized_keys\n\
if [ -n "$SSH_PUBLIC_KEY" ]; then\n\
  echo "$SSH_PUBLIC_KEY" > /root/.ssh/authorized_keys\n\
fi\n\
chmod 700 /root/.ssh\n\
chmod 600 /root/.ssh/authorized_keys\n\
exec /usr/sbin/sshd -D\n' > /etc/services.d/sshd/run && \
    chmod +x /etc/services.d/sshd/run

# 7. Expose the web port and SSH port
EXPOSE 8080 22

# (Note: No CMD is needed. The Rocker base image automatically starts the RStudio Server process)