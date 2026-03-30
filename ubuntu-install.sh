#!/bin/bash

# Methylation Analysis Environment Setup Script
# This script installs R, system dependencies, and sets up renv for the project

set -e  # Exit on error
set -o pipefail  # Exit if any command in a pipe fails

echo "=========================================="
echo "Setting up Methylation Analysis Environment"
echo "=========================================="

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Function to print colored output
print_status() {
    echo -e "${GREEN}[INFO]${NC} $1"
}

print_error() {
    echo -e "${RED}[ERROR]${NC} $1"
}

print_warning() {
    echo -e "${YELLOW}[WARNING]${NC} $1"
}

# Step 1: Update system and install R
print_status "Step 1: Installing R and system dependencies..."

# Add R repository for latest version
sudo apt-get update
sudo apt-get install -y software-properties-common dirmngr

# Add CRAN repository (optional but recommended for latest R)
# sudo add-apt-repository -y ppa:marutter/rrutter4
# sudo add-apt-repository -y ppa:c2d4u.team/c2d4u4.0+

# Update again after adding repositories
sudo apt-get update

# Install R and base packages
sudo apt-get install -y \
    r-base \
    r-base-dev \
    r-base-core \
    r-recommended

# # Step 2: Install system dependencies for Bioconductor packages
print_status "Step 2: Installing system dependencies for R packages..."

sudo apt-get install -y libcurl4-openssl-dev

sudo apt-get install -y \
    build-essential \
    gfortran \
    libssl-dev \
    libxml2-dev \
    libz-dev \
    libbz2-dev \
    liblzma-dev \
    libpcre2-dev \
    libgit2-dev \
    cmake \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    libfontconfig1-dev \
    libhdf5-dev \
    libfftw3-dev \
    libudunits2-dev \
    libgdal-dev \
    libgeos-dev \
    libproj-dev \
    libsqlite3-dev \
    libhts-dev \
    libncurses-dev \
    libblas-dev \
    liblapack-dev \
    libreadline-dev \
    libxt-dev \
    xorg-dev

if [ -f /usr/lib/x86_64-linux-gnu/libgfortran.so.5 ]; then
    sudo ln -sf /usr/lib/x86_64-linux-gnu/libgfortran.so.5 /usr/lib/x86_64-linux-gnu/libgfortran.so
    sudo ldconfig
fi

# Install libhts-dev separately with specific handling
print_status "Installing libhts-dev with openssl support..."

# Check if libhts-dev is already installed
if dpkg -l | grep -q libhts-dev; then
    print_status "libhts-dev is already installed"
else
    # Try to install libhts-dev allowing dependency resolution
    sudo apt-get install -y --no-install-recommends libhts-dev || \
    sudo apt-get install -y -f libhts-dev || \
    print_warning "libhts-dev installation had issues, but continuing..."
fi

# Verify libhts-dev installation
if dpkg -l | grep -q libhts-dev; then
    print_status "✓ libhts-dev installed successfully"
else
    print_warning "libhts-dev not properly installed. Rhtslib may fail to compile."
fi

# Step 3: Verify R installation
print_status "Step 3: Verifying R installation..."
R --version

if [ $? -ne 0 ]; then
    print_error "R installation failed!"
    exit 1
fi

# Step 4: Set up renv and install Bioconductor packages
print_status "Step 4: Setting up renv and installing Bioconductor packages..."

# Create R script for package installation
cat > setup_packages.R << 'EOF'
#!/usr/bin/env Rscript

# Set CRAN mirror
options(repos = c(CRAN = "https://cloud.r-project.org"))

# Print status
cat("Installing renv...\n")
install.packages("renv", quiet = FALSE)

cat("Initializing renv...\n")
renv::init(bare = TRUE)

# Set compiler flags for better compatibility
Sys.setenv(PKG_CFLAGS = "-O3")
Sys.setenv(PKG_CXXFLAGS = "-O3")
Sys.setenv(MAKEFLAGS = "-j2")  # Adjust based on your CPU cores

install.packages("progress", "tidyverse")
install.packages("patchwork")
install.packages("rlist")
install.packages("devtools")
install.packages("Rserve")

# Install BiocManager
cat("Installing BiocManager...\n")
install.packages("BiocManager", quiet = FALSE)

devtools::install_github("markgene/maxprobes")

# Install packages one by one to isolate failures
packages <- c(
    "Rhtslib",
    "Rhdf5lib",
    "IlluminaHumanMethylation450kanno.ilmn12.hg19",
    "IlluminaHumanMethylation450kmanifest",
    "minfi",
    "wateRmelon",
    "methylumi",
    "GEOquery",
    "limma",
    "ChAMP",
    "tidyverse",
    "data.table",
    "FlowSorted.Blood.450k",
    "FlowSorted.Blood.EPIC",
    "clusterProfiler",
    "pheatmap",
    "DMRcate",
    "reshape2",
    "RPMM",
    "ggpubr",
    "minfiData",
    "GenomicRanges",
    "BSgenome.Hsapiens.UCSC.hg19"
)

for (pkg in packages) {
    cat("\n========================================\n")
    cat("Installing:", pkg, "\n")
    cat("========================================\n")
    
    tryCatch({
        BiocManager::install(pkg, update = FALSE, ask = FALSE)
        cat("✓ Successfully installed:", pkg, "\n")
    }, error = function(e) {
        cat("✗ Failed to install:", pkg, "\n")
        cat("Error message:", e$message, "\n")
    })
}

# Take snapshot of installed packages
cat("\nCreating renv snapshot...\n")
renv::snapshot(confirm = FALSE)

cat("\n✓ Package installation completed!\n")

# List installed packages
cat("\nInstalled packages:\n")
installed <- renv::dependencies()
print(installed)

EOF

# Run the R setup script
print_status "Running R package installation script..."
R --vanilla < setup_packages.R

if [ $? -ne 0 ]; then
    print_error "Package installation failed!"
    exit 1
fi

# Step 5: Verify installation
print_status "Step 5: Verifying package installation..."

cat > verify_packages.R << 'EOF'
#!/usr/bin/env Rscript

# Function to check if package is installed
check_package <- function(pkg) {
    if (require(pkg, character.only = TRUE, quietly = TRUE)) {
        cat("✓", pkg, "- version:", as.character(packageVersion(pkg)), "\n")
        return(TRUE)
    } else {
        cat("✗", pkg, "- NOT INSTALLED\n")
        return(FALSE)
    }
}

# Check critical packages
critical_packages <- c(
    "minfi",
    "ChAMP",
    "GEOquery",
    "IlluminaHumanMethylation450kanno.ilmn12.hg19"
)

cat("\nVerifying critical packages:\n")
cat("========================================\n")

all_ok <- TRUE
for (pkg in critical_packages) {
    if (!check_package(pkg)) {
        all_ok <- FALSE
    }
}

if (all_ok) {
    cat("\n✓ All critical packages installed successfully!\n")
} else {
    cat("\n⚠ Some packages failed to install. Check the errors above.\n")
}

EOF

R --vanilla < verify_packages.R

# Step 6: Final setup instructions
print_status "Setup complete!"

cat << EOF

========================================
Setup Summary
========================================

✓ R installed: $(R --version | head -1)
✓ System dependencies installed
✓ renv initialized
✓ Bioconductor packages installed

Next steps:
1. Activate renv environment:
   R --vanilla -e "renv::activate()"

2. Run your Dependencies.R if needed:
   Rscript --vanilla Dependencies.R

3. Run your main analysis:
   Rscript --vanilla GSE145361.R

4. To deactivate renv (if needed):
   R --vanilla -e "renv::deactivate()"

Troubleshooting:
- If you encounter memory issues, try: export R_MAX_VSIZE=16Gb
- To see installed packages: R --vanilla -e "renv::dependencies()"
- To update packages: R --vanilla -e "renv::update()"

========================================
EOF

# Optional: Create a convenience script to run the analysis
cat > run_analysis.sh << 'EOF'
#!/bin/bash
# Convenience script to run the methylation analysis

echo "Activating renv and running GSE145361.R..."
Rscript --vanilla GSE145361.R

if [ $? -eq 0 ]; then
    echo "✓ Analysis completed successfully!"
else
    echo "✗ Analysis failed. Check the error messages above."
    exit 1
fi
EOF

chmod +x run_analysis.sh
print_status "Created convenience script: ./run_analysis.sh"

print_status "Environment setup complete! You can now run your analysis."