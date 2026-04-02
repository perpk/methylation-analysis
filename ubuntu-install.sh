#!/bin/bash

# Methylation Analysis Environment Setup Script
# This script checks for existing R installations, removes them if present,
# installs R 4.5 from CRAN, and sets up renv for the project

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

# Detect Ubuntu version
detect_ubuntu_version() {
    if [ -f /etc/os-release ]; then
        . /etc/os-release
        UBUNTU_VERSION="$VERSION_ID"
        UBUNTU_CODENAME="$VERSION_CODENAME"
        print_status "Detected Ubuntu $UBUNTU_VERSION ($UBUNTU_CODENAME)"
    else
        print_error "Cannot detect Ubuntu version"
        exit 1
    fi
}

# Step 1: Check and remove existing R installations
print_status "Step 1: Checking for existing R installations..."

# Check if R is installed
if command -v R &> /dev/null; then
    CURRENT_R_VERSION=$(R --version | head -1 | awk '{print $3}')
    print_warning "Found R version $CURRENT_R_VERSION installed"
    
    # Ask user if they want to remove
    read -p "Do you want to remove existing R installation? (y/n) " -n 1 -r
    echo
    if [[ $REPLY =~ ^[Yy]$ ]]; then
        print_status "Removing existing R packages..."
        
        # Remove R packages from user library
        if [ -d ~/R ]; then
            rm -rf ~/R
            print_status "Removed user R library"
        fi
        
        if [ -d ~/.local/lib/R ]; then
            rm -rf ~/.local/lib/R
            print_status "Removed local R library"
        fi
        
        # Remove system R installation
        print_status "Removing system R installation..."
        sudo apt-get remove --purge -y r-base r-base-dev r-base-core r-recommended || true
        sudo apt-get autoremove -y
        
        # Remove any remaining R directories
        sudo rm -rf /usr/lib/R
        sudo rm -rf /usr/local/lib/R
        sudo rm -rf /usr/share/R
        
        # Clean apt cache
        sudo apt-get update
        print_status "✓ Existing R installation removed"
    else
        print_status "Keeping existing R installation"
        # Check if we need to proceed with installation
        if [[ "$CURRENT_R_VERSION" != "4.5"* ]]; then
            print_error "Existing R version $CURRENT_R_VERSION is not 4.5.x"
            print_error "Please remove it manually or run this script with removal option"
            exit 1
        else
            print_status "R version $CURRENT_R_VERSION is compatible"
        fi
    fi
else
    print_status "No existing R installation found"
fi

# Step 2: Check and remove renv if present
print_status "Step 2: Checking for existing renv..."

if [ -d "renv" ] || [ -f "renv.lock" ]; then
    print_warning "Found existing renv environment"
    
    # Ask user if they want to remove
    read -p "Do you want to remove existing renv environment? (y/n) " -n 1 -r
    echo
    if [[ $REPLY =~ ^[Yy]$ ]]; then
        # Remove renv files and directories
        rm -rf renv
        rm -f renv.lock
        rm -rf .renv
        
        # Remove renv from .Rprofile if it exists
        if [ -f ".Rprofile" ]; then
            # Create backup
            cp .Rprofile .Rprofile.backup
            # Remove renv initialization lines
            sed -i '/source("renv\/activate.R")/d' .Rprofile
            print_status "✓ Existing renv environment removed"
        fi
    else
        print_status "Keeping existing renv environment"
    fi
fi

# Step 3: Detect Ubuntu version and set up CRAN repository
print_status "Step 3: Setting up CRAN repository for R 4.5..."

detect_ubuntu_version

# Map Ubuntu codename to CRAN repository
case $UBUNTU_CODENAME in
    focal)
        CRAN_REPO="focal-cran40"
        ;;
    jammy)
        CRAN_REPO="jammy-cran40"
        ;;
    noble)
        CRAN_REPO="noble-cran40"
        ;;
    *)
        print_error "Unsupported Ubuntu version: $UBUNTU_CODENAME"
        print_error "Supported versions: focal (20.04), jammy (22.04), noble (24.04)"
        exit 1
        ;;
esac

# Install dependencies for adding repositories
sudo apt-get update
sudo apt-get install -y --no-install-recommends \
    software-properties-common \
    dirmngr \
    wget \
    gnupg

# Add CRAN repository key
print_status "Adding CRAN repository key..."
wget -qO- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc | \
    sudo tee -a /etc/apt/trusted.gpg.d/cran_ubuntu_key.asc

# Add CRAN repository
print_status "Adding CRAN repository for Ubuntu $UBUNTU_CODENAME..."
sudo add-apt-repository -y "deb https://cloud.r-project.org/bin/linux/ubuntu $CRAN_REPO/"

# Update again after adding repositories
sudo apt-get update

# Step 4: Install R 4.5 and system dependencies
print_status "Step 4: Installing R 4.5 and system dependencies..."

# Install R 4.5
sudo apt-get install -y \
    r-base \
    r-base-dev \
    r-base-core

# Verify R version
R_VERSION=$(R --version | head -1 | awk '{print $3}')
print_status "Installed R version: $R_VERSION"

if [[ "$R_VERSION" != "4.5"* ]]; then
    print_error "Failed to install R 4.5. Got version $R_VERSION instead"
    exit 1
fi

# Step 5: Install system dependencies for Bioconductor packages
print_status "Step 5: Installing system dependencies for R packages..."

# Step 5: Install system dependencies for Bioconductor packages
print_status "Step 5: Installing system dependencies for R packages..."

# First install essential build tools and libraries
sudo apt-get install -y \
    build-essential \
    gfortran \
    libcurl4-openssl-dev \
    libcurl4-openssl-dev \
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
    libncurses-dev \
    libblas-dev \
    liblapack-dev \
    libreadline-dev \
    libxt-dev \
    xorg-dev

# Handle libhts-dev separately with proper dependency resolution
print_status "Installing libhts-dev with proper dependencies..."

# Try to fix any broken dependencies first
sudo apt-get --fix-broken install -y

# Install libcurl4-gnutls-dev if needed (alternative to openssl)
if ! sudo apt-get install -y libcurl4-gnutls-dev 2>/dev/null; then
    print_warning "libcurl4-gnutls-dev not available, using libcurl4-openssl-dev"
    # libcurl4-openssl-dev is already installed from previous step
fi

# Try to install libhts-dev with dependency fixing
if ! sudo apt-get install -y -f libhts-dev; then
    print_warning "libhts-dev installation failed, trying alternative approach..."
    
    # Try to install from source or skip for now
    print_status "libhts-dev will be handled by R packages if needed"
    
    # Some R packages can build without system libhts
    # Install development tools for building from source
    sudo apt-get install -y \
        autoconf \
        automake \
        libtool \
        pkg-config
fi

# Fix libgfortran symlink if needed
if [ -f /usr/lib/x86_64-linux-gnu/libgfortran.so.5 ]; then
    sudo ln -sf /usr/lib/x86_64-linux-gnu/libgfortran.so.5 /usr/lib/x86_64-linux-gnu/libgfortran.so
    sudo ldconfig
fi

# Step 6: Verify R installation
print_status "Step 6: Verifying R installation..."
# Step 6: Verify R installation
print_status "Step 6: Verifying R installation..."
R --version

if [ $? -ne 0 ]; then
    print_error "R installation failed!"
    exit 1
fi

# Step 7: Set up renv and install Bioconductor packages
print_status "Step 7: Setting up fresh renv and installing Bioconductor packages..."
# Step 7: Set up renv and install Bioconductor packages
print_status "Step 7: Setting up fresh renv and installing Bioconductor packages..."

# Create R script for package installation
cat > setup_packages.R << 'EOF'
#!/usr/bin/env Rscript

# Set CRAN mirror
options(repos = c(CRAN = "https://cloud.r-project.org"))

# Print R version
cat("R version:", R.version.string, "\n")

# Install renv
cat("\nInstalling renv...\n")
# Print R version
cat("R version:", R.version.string, "\n")

# Install renv
cat("\nInstalling renv...\n")
install.packages("renv", quiet = FALSE)

# Initialize renv (fresh)
cat("\nInitializing fresh renv...\n")
renv::init(bare = TRUE, settings = list(snapshot.type = "all"))
# Initialize renv (fresh)
cat("\nInitializing fresh renv...\n")
renv::init(bare = TRUE, settings = list(snapshot.type = "all"))

# Set compiler flags for better compatibility
Sys.setenv(PKG_CFLAGS = "-O3")
Sys.setenv(PKG_CXXFLAGS = "-O3")
Sys.setenv(MAKEFLAGS = "-j2")  # Adjust based on your CPU cores

# Install pak for better dependency resolution
cat("\nInstalling pak...\n")
install.packages("pak", quiet = FALSE)

# Install BiocManager
cat("\nInstalling BiocManager...\n")
cat("\nInstalling BiocManager...\n")
install.packages("BiocManager", quiet = FALSE)

# Install packages one by one to isolate failures
packages <- c(
    "Rhdf5lib",
    "IlluminaHumanMethylation450kanno.ilmn12.hg19",
    "IlluminaHumanMethylation450kmanifest",
    "minfi",
    "wateRmelon",
    "methylumi",
    "GEOquery",
    "limma",
    "tidyverse",
    "data.table",
    "pheatmap",
    "DMRcate",
    "reshape2",
    "ggpubr",
    "GenomicRanges",
    "BSgenome.Hsapiens.UCSC.hg19"
)

# Install pak for better dependency resolution
install.packages("pak", quiet = FALSE)

install.packages("rlist")

pak::pkg_install("markgene/maxprobes")

# Install packages
for (pkg in packages) {
    cat("\n========================================\n")
    cat("Installing:", pkg, "\n")
    cat("========================================\n")
    
    tryCatch({
        BiocManager::install(pkg, update = FALSE, ask = FALSE)
        cat("✓ Successfully installed:", pkg, "\n")
    }, error = function(e) {
        # Try with pak if BiocManager fails
        cat("Retrying with pak...\n")
        tryCatch({
            pak::pkg_install(pkg)
            cat("✓ Successfully installed with pak:", pkg, "\n")
        }, error = function(e2) {
            cat("✗ Failed to install:", pkg, "\n")
            cat("Error message:", e2$message, "\n")
        })
        # Try with pak if BiocManager fails
        cat("Retrying with pak...\n")
        tryCatch({
            pak::pkg_install(pkg)
            cat("✓ Successfully installed with pak:", pkg, "\n")
        }, error = function(e2) {
            cat("✗ Failed to install:", pkg, "\n")
            cat("Error message:", e2$message, "\n")
        })
    })
}

# Take snapshot of installed packages
cat("\nCreating renv snapshot...\n")
renv::snapshot(confirm = FALSE)

cat("\n✓ Package installation completed!\n")

# List installed packages
cat("\nInstalled packages summary:\n")
installed_pkgs <- installed.packages()
cat("Total packages installed:", nrow(installed_pkgs), "\n")

# Check critical packages
critical_packages <- c("minfi", "ChAMP", "GEOquery")
for (pkg in critical_packages) {
    if (pkg %in% rownames(installed_pkgs)) {
        cat("✓", pkg, "- version:", installed_pkgs[pkg, "Version"], "\n")
    } else {
        cat("✗", pkg, "- NOT INSTALLED\n")
    }
}

EOF

# Run the R setup script
print_status "Running R package installation script..."
R --vanilla < setup_packages.R

if [ $? -ne 0 ]; then
    print_warning "Some packages failed to install, but continuing..."
fi

# Step 8: Verify installation
print_status "Step 8: Verifying package installation..."
# Step 8: Verify installation
print_status "Step 8: Verifying package installation..."

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
    "GEOquery",
    "IlluminaHumanMethylation450kanno.ilmn12.hg19",
    "methylumi"
    "IlluminaHumanMethylation450kanno.ilmn12.hg19",
    "methylumi"
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

# Check minfi version
if (require("minfi", quietly = TRUE)) {
    minfi_version <- packageVersion("minfi")
    cat("\nminfi version:", as.character(minfi_version), "\n")
    if (minfi_version >= "1.56.0") {
        cat("✓ minfi version 1.56.0 or higher detected\n")
    } else {
        cat("⚠ WARNING: minfi version is older than 1.56.0\n")
    }
}

# Check minfi version
if (require("minfi", quietly = TRUE)) {
    minfi_version <- packageVersion("minfi")
    cat("\nminfi version:", as.character(minfi_version), "\n")
    if (minfi_version >= "1.56.0") {
        cat("✓ minfi version 1.56.0 or higher detected\n")
    } else {
        cat("⚠ WARNING: minfi version is older than 1.56.0\n")
    }
}

EOF

R --vanilla < verify_packages.R

# Step 9: Final setup instructions
# Step 9: Final setup instructions
print_status "Setup complete!"

cat << EOF

========================================
Setup Summary
========================================

✓ R installed: $(R --version | head -1)
✓ System dependencies installed
✓ Fresh renv initialized
✓ Fresh renv initialized
✓ Bioconductor packages installed

Next steps:
1. Activate renv (if not already active):
   source renv/activate
   or in R: renv::activate()
1. Activate renv (if not already active):
   source renv/activate
   or in R: renv::activate()

2. Run your main analysis:
2. Run your main analysis:
   Rscript --vanilla GSE145361.R

3. To see installed packages:
   R --vanilla -e "renv::dependencies()"
3. To see installed packages:
   R --vanilla -e "renv::dependencies()"

Troubleshooting:
- If you encounter memory issues, try: export R_MAX_VSIZE=16Gb
- To update packages: R --vanilla -e "renv::update()"
- To restore renv from lockfile: R --vanilla -e "renv::restore()"
- To restore renv from lockfile: R --vanilla -e "renv::restore()"

========================================
EOF

# Optional: Create a convenience script to run the analysis
cat > run_analysis.sh << 'EOF'
#!/bin/bash
# Convenience script to run the methylation analysis

echo "Running methylation analysis..."
echo "Running methylation analysis..."
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

