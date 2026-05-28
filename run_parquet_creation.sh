#!/bin/bash

PROJECT_NAME="GSE145361"
ROOT_DIR="/Volumes/Elements/vastai/combat/"
TARGET_DF_LOC="GSE145361_harmonized_targets.rds"
M_VALUES_LOC="GSE145361_combat_m_values.rds"
HARMONIZE_TARGETS="TRUE"
# for ppmi change to "ppmi/Project120_IDATS_n524final_toLONI_030718/"
# for GEO-datasets set to e.g. GSE111629_RAW/
PATTERN_HARM="GSE111629_RAW/"

TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
FILENAME="${PROJECT_NAME}_${TIMESTAMP}"
PID_FILE="${FILENAME}.pid"

nohup Rscript ./R/rds_to_parquet_converter.R "$PROJECT_NAME" "$ROOT_DIR" "$TARGET_DF_LOC" "$M_VALUES_LOC" "$HARMONIZE_TARGETS" "$PATTERN_HARM"  "$@" > "${FILENAME}.log" 2>&1 &

R_ID=$!

echo $R_ID > "$PID_FILE"