#!/usr/bin/env python3
"""
Amalgkit Status Monitor (MCP Tool).

This script functions as a Model Context Protocol (MCP) tool.
It parses the Amalgkit pipeline logs and process table to return 
structured JSON status information for AI agents.

Usage:
    python3 -m metainformant.mcp.tools.amalgkit_monitor
    
Output (JSON):
    {
        "status": "running" | "stalled" | "completed" | "error",
        "pid": 12345,
        "workers": 20,
        "progress": {
            "processed": 35,
            "total": 5615,
            "percent": 0.62
        },
        "last_log": "[35/5615] ✓ SRR26150181 (133.1s)...",
        "system": {
            "load_avg": [15.5, 12.0, 10.0],
            "free_disk_gb": 450.2
        }
    }
"""

import json
import os
import re
import sys
import psutil
from pathlib import Path

# Configuration - could be dynamic in future
WORK_DIR = Path("/Volumes/blue/data/amalgkit/apis_mellifera_all/work")
LOG_FILE = WORK_DIR / "ena_parallel_processing.log"

def get_process_status():
    """Find the active Amalgkit process."""
    for proc in psutil.process_iter(['pid', 'name', 'cmdline']):
        try:
            cmdline = proc.info['cmdline'] or []
            # Check for our specific script name
            if any("process_apis_mellifera" in arg for arg in cmdline) or \
               any("process_species.py" in arg for arg in cmdline):
                return {
                    "running": True,
                    "pid": proc.info['pid'],
                    "cmd": " ".join(cmdline)
                }
        except (psutil.NoSuchProcess, psutil.AccessDenied):
            continue
    return {"running": False, "pid": None}

def parse_log_progress():
    """Parse the last relevant line from the log file."""
    if not LOG_FILE.exists():
        return {"processed": 0, "total": 0, "last_line": "No log file found"}

    try:
        # Read last few lines to find progress
        # Using tail-like logic
        with open(LOG_FILE, 'rb') as f:
            try:
                f.seek(-2000, os.SEEK_END)
            except OSError:
                f.seek(0)
            
            lines = f.readlines()
            
            # Decode and look for progress pattern: [35/5615]
            # Pattern: matches [123/456] key
            progress_re = re.compile(r'\[(\d+)/(\d+)\]')
            
            processed = 0
            total = 0
            last_line = ""

            for line in reversed(lines):
                line_str = line.decode('utf-8', errors='ignore').strip()
                if not line_str:
                    continue
                
                match = progress_re.search(line_str)
                if match:
                    processed = int(match.group(1))
                    total = int(match.group(2))
                    last_line = line_str
                    break
            
            return {
                "processed": processed,
                "total": total,
                "percent": round((processed / total) * 100, 2) if total > 0 else 0,
                "last_line": last_line
            }

    except Exception as e:
        return {"error": str(e)}

def get_system_stats():
    """Get system load and disk space."""
    load = os.getloadavg()
    
    # Check disk space on work dir volume
    try:
        stat = os.statvfs(WORK_DIR)
        free_gb = (stat.f_frsize * stat.f_bavail) / (1024**3)
    except OSError:
        free_gb = 0

    return {
        "load_avg": load,
        "free_disk_gb": round(free_gb, 1)
    }

def main():
    """Main execution."""
    proc_info = get_process_status()
    log_info = parse_log_progress()
    sys_info = get_system_stats()

    status = "stopped"
    if proc_info["running"]:
        status = "running"
        # Heuristic for stalled: running but load is very low? 
        # Or simple check: is it running?
    
    output = {
        "status": status,
        "process": proc_info,
        "progress": {
            "processed": log_info.get("processed", 0),
            "total": log_info.get("total", 0),
            "percent": log_info.get("percent", 0)
        },
        "last_log": log_info.get("last_line", ""),
        "system": sys_info
    }

    print(json.dumps(output, indent=2))

if __name__ == "__main__":
    main()
