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
install.packages("patchwork")

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
    "BSgenome.Hsapiens.UCSC.hg19",
    "FlowSorted.Blood.450k",
    "FlowSorted.Blood.EPIC",
    "RPMM",
    "sva"
)

# Install pak for better dependency resolution
install.packages("pak", quiet = FALSE)

install.packages("rlist")
install.packages("arrow")

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

