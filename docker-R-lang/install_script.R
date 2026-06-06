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
    "missMethyl",
    "clusterProfiler",
    "sva",
    "minfiData",
    "annotatr",
    "GenomicRanges",
    "BSgenome.Hsapiens.UCSC.hg19"
)

install.packages("progress", "tidyverse", "BiocManager")
# You need that for the BMIQ plot.
install.packages("patchwork")
install.packages("devtools")
install.packages("Rserve")

# Install pak for better dependency resolution
install.packages("pak", quiet = FALSE)

install.packages("vscDebugger", repos = "https://manuelhentschel.r-universe.dev")
devtools::install_github("markgene/maxprobes")
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
