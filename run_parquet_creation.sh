#!/bin/bash

PROJECT_NAME="GSE111629"
ROOT_DIR="/root/workspace/methyl-pipe-out/GSE111629_20260515_083253/"
TARGET_DF_LOC="processed/targets_remove_mismatch.rds"
M_VALUES_LOC="results/m_values_bmiq.rds"
HARMONIZE_TARGETS="TRUE"
# for ppmi change to ^\\d{4}_
PATTERN_HARM="GSE111629_RAW/"

TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
FILENAME="${PROJECT_NAME}_${TIMESTAMP}"
PID_FILE="${FILENAME}.pid"

nohup Rscript ./R/rds_to_parquet_converter.R "$PROJECT_NAME" "$ROOT_DIR" "$TARGET_DF_LOC" "$M_VALUES_LOC" "$HARMONIZE_TARGETS" "$PATTERN_HARM"  "$@" > "${FILENAME}.log" 2>&1 &

R_ID=$!

echo $R_ID > "$PID_FILE"