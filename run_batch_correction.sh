#!/bin/bash

PROJECT_NAME <- ""
ROOT_DIR <- ""
TARGET_DF_LOC <- ""
M_VALUES_LOC <- ""
IDAT_FOLDER_LOC <- ""

TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
FILENAME="${PROJECT_NAME}_${TIMESTAMP}"
PID_FILE="${PROJECT_NAME}.pid"
nohup Rscript ./R/batch_correction_central.R "$PROJECT_NAME" "$ROOT_DIR" "$TARGET_DF_LOC" "$M_VALUES_LOC" "$IDAT_FOLDER_LOC" "$@" > "${FILENAME}.log" 2>&1 &

R_PID=$!

echo $R_PID > "$PID_FILE"