#!/bin/bash
# Monitor script for GSE111629.R (PID: 400236)

echo "Monitoring R process (PID: 400236)"
echo "Log file: logs/GSE111629_20260331_141633.log"
echo ""
echo "Recent log entries:"
tail -20 "logs/GSE111629_20260331_141633.log"
echo ""
echo "Process status:"
ps -p 400236 -o pid,ppid,state,etime,cmd 2>/dev/null || echo "Process not running"
echo ""
echo "To follow logs: tail -f logs/GSE111629_20260331_141633.log"
echo "To kill process: kill 400236"
