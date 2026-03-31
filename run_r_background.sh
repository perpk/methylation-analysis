#!/bin/bash

# run_r_background.sh - Run an R script in the background with nohup and logging
# Usage: ./run_r_background.sh <R_script> [additional Rscript arguments]

set -e  # Exit on error

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Function to print colored output
print_status() {
    echo -e "${GREEN}[INFO]${NC} $1"
}

print_error() {
    echo -e "${RED}[ERROR]${NC} $1"
}

print_warning() {
    echo -e "${YELLOW}[WARNING]${NC} $1"
}

# Check if R script is provided
if [ $# -eq 0 ]; then
    print_error "No R script specified"
    echo "Usage: $0 <R_script.R> [Rscript arguments...]"
    echo ""
    echo "Examples:"
    echo "  $0 analysis.R"
    echo "  $0 analysis.R --args --option value"
    echo "  $0 /path/to/script.R param1 param2"
    exit 1
fi

R_SCRIPT="$1"
shift  # Remove the first argument, leaving any additional arguments

# Check if R script exists
if [ ! -f "$R_SCRIPT" ]; then
    print_error "R script not found: $R_SCRIPT"
    exit 1
fi

# Check if Rscript is available
if ! command -v Rscript &> /dev/null; then
    print_error "Rscript not found. Please install R first."
    exit 1
fi

# Create logs directory if it doesn't exist
LOGS_DIR="logs"
mkdir -p "$LOGS_DIR"

# Generate timestamp for log file
TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
SCRIPT_NAME=$(basename "$R_SCRIPT" .R)
LOG_FILE="${LOGS_DIR}/${SCRIPT_NAME}_${TIMESTAMP}.log"
PID_FILE="${LOGS_DIR}/${SCRIPT_NAME}_${TIMESTAMP}.pid"

print_status "R Script: $R_SCRIPT"
print_status "Log file: $LOG_FILE"
print_status "PID file: $PID_FILE"
print_status "Arguments: $@"

# Run R script with nohup in background
print_status "Starting R script in background..."

nohup Rscript "$R_SCRIPT" "$@" > "$LOG_FILE" 2>&1 &

# Capture the process ID
R_PID=$!

# Save PID to file
echo $R_PID > "$PID_FILE"

print_status "R process started with PID: $R_PID"
print_status "PID saved to: $PID_FILE"
print_status "Log file: $LOG_FILE"

# Optional: Wait a moment and check if process is still running
sleep 2

if kill -0 $R_PID 2>/dev/null; then
    print_status "✓ Process is running successfully"
    print_status ""
    print_status "Monitor progress with:"
    print_status "  tail -f $LOG_FILE"
    print_status ""
    print_status "Check if process is still running:"
    print_status "  ps -p $R_PID"
    print_status ""
    print_status "Kill the process if needed:"
    print_status "  kill $R_PID"
else
    print_warning "Process may have failed to start. Check log file:"
    print_warning "  cat $LOG_FILE"
    exit 1
fi

# Create a helper script to monitor this specific job
MONITOR_SCRIPT="${LOGS_DIR}/monitor_${SCRIPT_NAME}_${TIMESTAMP}.sh"
cat > "$MONITOR_SCRIPT" << EOF
#!/bin/bash
# Monitor script for $R_SCRIPT (PID: $R_PID)

echo "Monitoring R process (PID: $R_PID)"
echo "Log file: $LOG_FILE"
echo ""
echo "Recent log entries:"
tail -20 "$LOG_FILE"
echo ""
echo "Process status:"
ps -p $R_PID -o pid,ppid,state,etime,cmd 2>/dev/null || echo "Process not running"
echo ""
echo "To follow logs: tail -f $LOG_FILE"
echo "To kill process: kill $R_PID"
EOF

chmod +x "$MONITOR_SCRIPT"
print_status "Created monitor script: $MONITOR_SCRIPT"

print_status ""
print_status "Setup complete! The R script is running in the background."
