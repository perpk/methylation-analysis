#!/bin/bash
# Convenience script to run the methylation analysis

echo "Running methylation analysis..."
Rscript --vanilla GSE145361.R

if [ $? -eq 0 ]; then
    echo "✓ Analysis completed successfully!"
else
    echo "✗ Analysis failed. Check the error messages above."
    exit 1
fi
