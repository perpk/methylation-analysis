#!/bin/bash

PROJECT_NAME="PPMI"
ROOT_DIR="/root/workspace/methyl-pipe-out//ppmi_20260513_110353/"
TARGET_DF_LOC="processed/targets_remove_mismatch.rds"
M_VALUES_LOC="results/m_values_bmiq.rds"
IDAT_FOLDER_LOC="/root/methylation-analysis/ppmi/Project120_IDATS_n524final_toLONI_030718"
EXTRACT_SENTRIX_ID_FROM_BASENAME="TRUE"
HARMONIZE_TARGETS="TRUE"
# "GSM\\d+_\\d{10}_R\\d{2}C\\d{2}"
PREFIX_PTN="\\d{12}_R\\d{2}C\\d{2}"

TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
FILENAME="${PROJECT_NAME}_${TIMESTAMP}"
PID_FILE="${FILENAME}.pid"
nohup Rscript ./R/batch_correction_central.R "$PROJECT_NAME" "$ROOT_DIR" "$TARGET_DF_LOC" "$M_VALUES_LOC" "$IDAT_FOLDER_LOC" "$EXTRACT_SENTRIX_ID_FROM_BASENAME" "$HARMONIZE_TARGETS" "$PREFIX_PTN" "$@" > "${FILENAME}.log" 2>&1 &

R_PID=$!

echo $R_PID > "$PID_FILE"