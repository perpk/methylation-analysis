#!/bin/bash
# Monitor script for GSE111629.R (PID: 331309)

echo "Monitoring R process (PID: 331309)"
echo "Log file: logs/GSE111629_20260331_073545.log"
echo ""
echo "Recent log entries:"
tail -20 "logs/GSE111629_20260331_073545.log"
echo ""
echo "Process status:"
ps -p 331309 -o pid,ppid,state,etime,cmd 2>/dev/null || echo "Process not running"
echo ""
echo "To follow logs: tail -f logs/GSE111629_20260331_073545.log"
echo "To kill process: kill 331309"
