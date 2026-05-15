#!/bin/bash

PROJECT_NAME="GSE111629"
ROOT_DIR="/root/workspace/methyl-pipe-out/GSE111629_20260515_083253/"
TARGET_DF_LOC="processed/targets_remove_mismatch.rds"
M_VALUES_LOC="results/m_values_bmiq.rds"
IDAT_FOLDER_LOC="/root/methylation-analysis/GSE111629_RAW"
EXTRACT_SENTRIX_ID_FROM_BASENAME="TRUE"
HARMONIZE_TARGETS="TRUE"

TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
FILENAME="${PROJECT_NAME}_${TIMESTAMP}"
PID_FILE="${FILENAME}.pid"
nohup Rscript ./R/batch_correction_central.R "$PROJECT_NAME" "$ROOT_DIR" "$TARGET_DF_LOC" "$M_VALUES_LOC" "$IDAT_FOLDER_LOC" "$EXTRACT_SENTRIX_ID_FROM_BASENAME" "$HARMONIZE_TARGETS" "$@" > "${FILENAME}.log" 2>&1 &

R_PID=$!

echo $R_PID > "$PID_FILE"