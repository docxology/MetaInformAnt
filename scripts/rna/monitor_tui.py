#!/usr/bin/env python3
"""
Amalgkit Pipeline TUI Monitor
=============================

Real-time visualization of the Amalgkit pipeline progress w/ System Metrics.
Shows:
- System: CPU, RAM, Network I/O, Active Commands
- Species: Active (Running), Downloaded, Quantified, Total
- Status: Merge, Curate, Sanity steps

Usage:
    python3 scripts/rna/monitor_tui.py
"""

import time
import os
import sys
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime

import psutil
from rich.live import Live
from rich.table import Table
from rich.console import Console, Group
from rich.panel import Panel
from rich.layout import Layout
from rich import box

# Constants
SPECIES_COUNTS = {
    "anoplolepis_gracilipes": 2, "acromyrmex_echinatior": 8, "dinoponera_quadriceps": 13,
    "vollenhovia_emeryi": 15, "odontomachus_brunneus": 19, "formica_exsecta": 23,
    "temnothorax_americanus": 32, "wasmannia_auropunctata": 33, "nylanderia_fulva": 40,
    "temnothorax_curvispinosus": 43, "pbarbatus": 95, "cardiocondyla_obscurior": 162,
    "temnothorax_nylanderi": 166, "linepithema_humile": 173, "atta_cephalotes": 220,
    "ooceraea_biroi": 237, "camponotus_floridanus": 304, "solenopsis_invicta": 349,
    "monomorium_pharaonis": 370, "temnothorax_longispinosus": 508,
    "harpegnathos_saltator": 689, "apis_mellifera": 3154, "amellifera": 3154
}

SPECIES_ORDER = [
    "anoplolepis_gracilipes", "acromyrmex_echinatior", "dinoponera_quadriceps",
    "vollenhovia_emeryi", "odontomachus_brunneus", "formica_exsecta",
    "temnothorax_americanus", "wasmannia_auropunctata", "nylanderia_fulva",
    "temnothorax_curvispinosus", "pbarbatus", "cardiocondyla_obscurior",
    "temnothorax_nylanderi", "linepithema_humile", "atta_cephalotes",
    "ooceraea_biroi", "camponotus_floridanus", "solenopsis_invicta",
    "monomorium_pharaonis", "temnothorax_longispinosus", "harpegnathos_saltator",
    "apis_mellifera"
]

AMALGKIT_DIR = Path("output/amalgkit")

class SystemMonitor:
    def __init__(self):
        self.last_net = psutil.net_io_counters()
        self.last_time = time.time()
        self.all_relevant_procs = []
        
    def get_metrics(self):
        # 1. Network Speed
        current_net = psutil.net_io_counters()
        current_time = time.time()
        elapsed = current_time - self.last_time
        
        if elapsed > 0:
            recv_speed = (current_net.bytes_recv - self.last_net.bytes_recv) / elapsed
            sent_speed = (current_net.bytes_sent - self.last_net.bytes_sent) / elapsed
        else:
            recv_speed = 0
            sent_speed = 0
            
        self.last_net = current_net
        self.last_time = current_time
        
        # 2. CPU / RAM (Total System)
        cpu_pct = psutil.cpu_percent(interval=None)
        mem = psutil.virtual_memory()
        
        # 3. Active Processes (Scan once)
        relevant_procs = []
        self.all_relevant_procs = [] # Clear cache
        
        target_names = {"amalgkit", "python", "python3", "kallisto", "prefetch", "fastq-dump", "curl", "wget", "fasterq-dump"}
        
        # Iterate over processes safely
        for proc in psutil.process_iter(['pid', 'name', 'cpu_percent', 'memory_info', 'cmdline']):
            try:
                cmdline = proc.info['cmdline'] or []
                cmd_str = " ".join(cmdline)
                
                if proc.info['name'] in target_names or any(t in cmd_str for t in target_names):
                    if "monitor_tui.py" in cmd_str:
                        continue
                        
                    mem_info = proc.info.get('memory_info')
                    mem_usage = 0.0
                    if mem_info:
                        mem_usage = mem_info.rss / (1024*1024) # MB

                    p_info = {
                        "pid": proc.info['pid'],
                        "name": proc.info['name'],
                        "cpu": proc.info['cpu_percent'],
                        "mem": mem_usage,
                        "cmd": cmd_str
                    }
                    self.all_relevant_procs.append(p_info)
            except (psutil.NoSuchProcess, psutil.AccessDenied, psutil.ZombieProcess, OSError):
                pass
                
        # Sort by CPU usage for display
        relevant_procs.sort(key=lambda x: x['cpu'], reverse=True)
        
        return {
            "net_recv": recv_speed,
            "net_sent": sent_speed,
            "cpu_total": cpu_pct,
            "mem_pct": mem.percent,
            "mem_used": mem.used / (1024**3), # GB
            "procs": relevant_procs[:5] # Top 5
        }
        
    def count_active_for_species(self, species_name):
        # Check cmdline for species name
        count = 0
        for p in self.all_relevant_procs:
            if species_name in p['cmd']:
                count += 1
        return count

def format_speed(bytes_per_sec):
    if bytes_per_sec > 1024**2:
        return f"{bytes_per_sec / (1024**2):.1f} MB/s"
    elif bytes_per_sec > 1024:
        return f"{bytes_per_sec / 1024:.1f} KB/s"
    else:
        return f"{bytes_per_sec:.0f} B/s"

def get_species_status(args):
    """Worker function needs to unpack args if using map."""
    species_name, active_count = args
    
    # Paths
    sp_dir = AMALGKIT_DIR / species_name
    if not sp_dir.exists():
        # Fallback to local output for logs/metadata if not in blue
        # But data should be in blue for 160GB download
        # TUI needs to check BOTH locations
        pass
        
    local_log_dir = Path("output/amalgkit") / species_name
    
    base_status = {
        "name": species_name,
        "active": active_count,
        "quant": 0,
        "downloaded": 0,
        "total": SPECIES_COUNTS.get(species_name, 100),
        "status": "Waiting",
        "merge": "[dim]·[/dim]",
        "curate": "[dim]·[/dim]",
        "sanity": "[dim]·[/dim]",
    }

    if not sp_dir.exists():
        return base_status

    # 1. Downloaded Count (Fastq files)
    # Location: blue/amalgkit/{species}/fastq/getfastq/{SRR}/{SRR}.fastq.gz
    # OR blue/amalgkit/{species}/work/getfastq/
    # Amalgkit default is usually defined in config.
    # Usually: fastq/getfastq per spec
    # Let's check typical paths
    fastq_dirs = [
        sp_dir / "fastq" / "getfastq",
        sp_dir / "work" / "getfastq"
    ]
    
    dl_count = 0
    for fd in fastq_dirs:
        if fd.exists():
            try:
                # Count directories that contain fastq.gz
                # Just counting directories in getfastq is a good proxy for 'attempted/downloaded'
                items = [x for x in fd.iterdir() if x.is_dir()]
                dl_count = len(items)
                if dl_count > 0: break
            except OSError:
                pass

    # 2. Quant Count
    quant_count = 0
    work_dir = sp_dir / "work"
    quant_dir = work_dir / "quant"
    if quant_dir.exists():
        try:
             items = [x for x in quant_dir.iterdir() if x.is_dir()]
             quant_count = len(items)
        except OSError:
            quant_count = 0
            
    # Refine Downloaded: Downloaded usually >= Quantified.
    # If dl_count < quant_count (e.g. clean_fastq=yes), then dl_count is ephemeral.
    # But usually valid.
    
    # 3. Flags
    has_merge = (sp_dir / "merged" / "merged_abundance.tsv").exists()
    curate_dir = work_dir / "curate"
    has_curate = False
    if curate_dir.exists():
        if list(curate_dir.glob("*_curated.tsv")) or (curate_dir / "curate_completion_flag.txt").exists():
            has_curate = True
    has_sanity = (work_dir / "sanity_completion_flag.txt").exists()
    
    # Status text
    status = "Waiting"
    if has_sanity:
        status = "[green]Done[/green]"
    elif has_curate:
        status = "[cyan]Curating[/cyan]"
    elif has_merge:
        status = "[blue]Merging[/blue]"
    elif active_count > 0:
        status = "[yellow]Running[/yellow]"
    elif quant_count > 0:
        # If not active but has quants, maybe paused or between steps?
        # Or fully done but flags missing?
        if quant_count >= base_status['total']:
             status = "[blue]Finishing[/blue]"
        else:
             status = "Paused"
    elif quant_count == 0 and sp_dir.exists():
        status = "Initializing"

    return {
        "name": species_name,
        "active": active_count,
        "quant": quant_count,
        "downloaded": dl_count,
        "total": SPECIES_COUNTS.get(species_name, 100),
        "status": status,
        "merge": "[green]Done[/green]" if has_merge else "[dim]·[/dim]",
        "curate": "[green]Done[/green]" if has_curate else "[dim]·[/dim]",
        "sanity": "[green]Done[/green]" if has_sanity else "[dim]·[/dim]",
    }

def generate_dashboard(sys_monitor, species_data) -> Group:
    metrics = sys_monitor.get_metrics()
    
    # 1. System Panel
    sys_table = Table.grid(expand=True, padding=(0, 2))
    sys_table.add_column("Metric", style="bold cyan")
    sys_table.add_column("Value", style="yellow")
    sys_table.add_row("CPU Usage", f"{metrics['cpu_total']}%")
    sys_table.add_row("Memory", f"{metrics['mem_used']:.1f} GB ({metrics['mem_pct']}%)")
    sys_table.add_row("Network (In)", format_speed(metrics['net_recv']))
    sys_table.add_row("Network (Out)", format_speed(metrics['net_sent']))
    
    proc_table = Table(box=box.SIMPLE_HEAD, title="Top Active Processes", title_style="bold magenta")
    proc_table.add_column("PID", style="dim")
    proc_table.add_column("Command", style="green")
    proc_table.add_column("CPU%", justify="right")
    
    for proc in metrics['procs']:
        proc_table.add_row(
            str(proc['pid']), 
            proc['cmd'][:50], 
            f"{proc['cpu']:.1f}"
        )
        
    sys_panel = Panel(
        Group(sys_table, proc_table),
        title=f"System Monitor - {datetime.now().strftime('%H:%M:%S')}",
        border_style="blue"
    )

    # 2. Species Table
    sp_table = Table(box=box.SIMPLE)
    sp_table.add_column("Species", style="cyan")
    sp_table.add_column("Status", style="bold")
    sp_table.add_column("Active", justify="right", style="yellow")
    sp_table.add_column("Dl/Quant/Total", justify="right")
    sp_table.add_column("Progress", style="magenta")
    sp_table.add_column("M/C/S", justify="center")

    for res in species_data:
        total = res['total']
        quant = res['quant']
        dl = res['downloaded']
        pct = (quant / total) * 100 if total > 0 else 0
        
        # Color bar
        bar_char = "█"
        width = 8
        filled = int((pct / 100) * width)
        bar = f"[{'green' if pct >= 100 else 'yellow'}]{bar_char * filled}[/][dim]{bar_char * (width - filled)}[/]"
        
        steps_str = f"{res['merge']}{res['curate']}{res['sanity']}".replace("Done", "✓")
        
        sp_table.add_row(
            res['name'],
            res['status'],
            str(res['active']) if res['active'] > 0 else "-",
            f"{dl}/{quant}/{total}",
            f"{bar} {pct:.0f}%",
            steps_str
        )

    sp_panel = Panel(sp_table, title="Pipeline Progress", border_style="green")
    
    return Group(sys_panel, sp_panel)

def main():
    console = Console()
    console.clear()
    
    monitor = SystemMonitor()
    time.sleep(0.1) # Init psutil
    
    with Live(refresh_per_second=1, console=console) as live:
        try:
            while True:
                # 1. Update system metrics & snapshot processes
                metrics = monitor.get_metrics() # This updates self.all_relevant_procs
                
                # 2. Prepare args for species status (name, active_count)
                status_args = []
                for sp in SPECIES_ORDER:
                    active = monitor.count_active_for_species(sp)
                    status_args.append((sp, active))
                
                # 3. Fetch status in thread pool
                with ThreadPoolExecutor(max_workers=8) as executor:
                    species_data = list(executor.map(get_species_status, status_args))
                
                dashboard = generate_dashboard(monitor, species_data)
                live.update(dashboard)
                time.sleep(2)
        except KeyboardInterrupt:
            console.print("\n[yellow]Monitor exited.[/yellow]")

if __name__ == "__main__":
    main()
