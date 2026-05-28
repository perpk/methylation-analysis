#!/bin/bash

PROJECTS="PPMI"
ROOT_DIR="/root/workspace/"

TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
FILENAME="BATCH_CORRECTION_${TIMESTAMP}"
PID_FILE="${FILENAME}.pid"
nohup Rscript ./R/batch_correction_central.R "$PROJECTS" "$ROOT_DIR" "$@" > "${FILENAME}.log" 2>&1 &

R_PID=$!

echo $R_PID > "$PID_FILE"
