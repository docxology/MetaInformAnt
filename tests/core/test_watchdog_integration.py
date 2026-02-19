"""
Integration test for ProcessWatchdog.
"""

import time
import subprocess
import sys
import pytest
from pathlib import Path
from metainformant.core.utils.watchdog import ProcessWatchdog

def test_watchdog_kills_stalled_process():
    """
    Test that watchdog kills a process that sleeps (low CPU) for too long.
    """
    # Create a python script that sleeps forever
    script = """
import time
import sys
print("Starting sleep...", flush=True)
try:
    while True:
        time.sleep(1)
except BaseException as e:
    print(f"Caught {type(e).__name__}", flush=True)
"""
    script_path = Path("temp_sleep_script.py")
    script_path.write_text(script)

    try:
        # Start the sleeper
        proc = subprocess.Popen(
            [sys.executable, str(script_path)],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True
        )
        
        # Start watchdog with very short timeout
        watchdog = ProcessWatchdog(
            pid=proc.pid,
            cpu_threshold=1.0, # Sleep uses ~0% CPU
            timeout_seconds=2, # Kill after 2 seconds
            check_interval=1,
            on_stall=ProcessWatchdog.kill_process_tree
        )
        watchdog.start()
        
        # Wait for potential kill
        try:
            # Should be killed within ~3-4 seconds
            stdout, stderr = proc.communicate(timeout=10)
            
            # If we get here, process ended. Check return code.
            # Killed processes usually have -9 (SIGKILL) or -15 (SIGTERM)
            assert proc.returncode != 0, "Process should have been killed (non-zero exit)"
            print(f"Process killed with rc={proc.returncode}")
            
        except subprocess.TimeoutExpired:
            proc.kill()
            pytest.fail("Watchdog failed to kill process within timeout")
            
        finally:
            watchdog.stop()

    finally:
        if script_path.exists():
            script_path.unlink()

if __name__ == "__main__":
    # Manual run support
    try:
        test_watchdog_kills_stalled_process()
        print("Test PASSED")
    except Exception as e:
        print(f"Test FAILED: {e}")
        sys.exit(1)
