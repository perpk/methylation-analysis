#!/bin/bash

PROJECT_NAME="PPMI"
ROOT_DIR="/root/workspace/methyl-pipe-out/ppmi_20260513_110353/"
TARGET_DF_LOC="processed/targets_remove_mismatch.rds"
M_VALUES_LOC="PPMI_combat_m_values.rds"
HARMONIZE_TARGETS="TRUE"
# for ppmi change to "ppmi/Project120_IDATS_n524final_toLONI_030718/200973410159_R03C01"
# for GEO-datasets set to e.g. GSE111629_RAW/
PATTERN_HARM="ppmi/Project120_IDATS_n524final_toLONI_030718/200973410159_R03C01"

TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
FILENAME="${PROJECT_NAME}_${TIMESTAMP}"
PID_FILE="${FILENAME}.pid"

nohup Rscript ./R/rds_to_parquet_converter.R "$PROJECT_NAME" "$ROOT_DIR" "$TARGET_DF_LOC" "$M_VALUES_LOC" "$HARMONIZE_TARGETS" "$PATTERN_HARM"  "$@" > "${FILENAME}.log" 2>&1 &

R_ID=$!

echo $R_ID > "$PID_FILE"