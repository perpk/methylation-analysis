# 1. Define packages (Duplicates removed, ChAMP added)
packages <- c(
    "Rhdf5lib",
    "IlluminaHumanMethylation450kanno.ilmn12.hg19",
    "IlluminaHumanMethylation450kmanifest",
    "minfi",
    "wateRmelon",
    "methylumi",
    "GEOquery",
    "limma",
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
    "ChAMP"
)

# 2. Base CRAN Installations
# Wrapped in c() and provided a default mirror for non-interactive execution
install.packages(
    c("progress", "BiocManager", "patchwork", "devtools", "Rserve", "pak", "rlist", "arrow", "renv", "mirai", "tidyverse"),
    repos = "http://cran.us.r-project.org",
    quiet = FALSE
)

# Custom repo installation
install.packages("vscDebugger", repos = "https://manuelhentschel.r-universe.dev")

# 3. GitHub Installations
pak::pkg_install("markgene/maxprobes", ask = FALSE)

# 4. Bulk Install Bioconductor/CRAN packages (Vectorized for speed)
cat("\n========================================\n")
cat("Starting bulk installation with BiocManager...\n")
cat("========================================\n")
BiocManager::install(packages, update = FALSE, ask = FALSE)

# 5. Smart Fallback using pak
# Check what is currently installed vs what was requested
installed_now <- rownames(installed.packages())
missing_pkgs <- setdiff(packages, installed_now)

if (length(missing_pkgs) > 0) {
    cat("\n========================================\n")
    cat("Attempting to install missing packages with pak:\n")
    print(missing_pkgs)
    cat("========================================\n")
    tryCatch({
        pak::pkg_install(missing_pkgs, ask = FALSE)
    }, error = function(e) {
        cat("✗ Some packages still failed to install via pak.\n")
        cat("Error message:", e$message, "\n")
    })
}

# 6. Take snapshot of installed packages
cat("\nCreating renv snapshot...\n")
# Must initialize the renv project first before snapping
if (!file.exists("renv.lock")) {
  renv::init(bare = TRUE, restart = FALSE)
}
renv::snapshot(confirm = FALSE)

cat("\n✓ Package installation completed!\n")

# 7. List installed packages summary
cat("\nInstalled packages summary:\n")
final_installed <- installed.packages()
cat("Total packages installed:", nrow(final_installed), "\n")

# 8. Check critical packages
critical_packages <- c("minfi", "ChAMP", "GEOquery")
for (pkg in critical_packages) {
    if (pkg %in% rownames(final_installed)) {
        cat("✓", pkg, "- version:", final_installed[pkg, "Version"], "\n")
    } else {
        cat("✗", pkg, "- NOT INSTALLED\n")
    }
}