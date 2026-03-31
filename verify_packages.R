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

