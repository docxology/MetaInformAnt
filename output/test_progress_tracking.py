#!/usr/bin/env python3
"""Test script for download progress tracking functionality."""

import sys
import tempfile
import time
from pathlib import Path

# Add src to path
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from metainformant.rna.steps.download_progress import (
    DownloadProgressMonitor,
    ThreadProgressTracker,
    FileSizeMonitor,
)


def test_file_size_monitor():
    """Test FileSizeMonitor functionality."""
    print("Testing FileSizeMonitor...")
    
    with tempfile.TemporaryDirectory() as tmpdir:
        test_dir = Path(tmpdir) / "test_sample"
        test_dir.mkdir()
        
        monitor = FileSizeMonitor(test_dir, update_interval=0.5)
        
        # Initial size should be 0
        assert monitor.get_size() == 0, "Initial size should be 0"
        print("  ✓ Initial size check passed")
        
        # Create a test file
        test_file = test_dir / "test.fastq.gz"
        test_file.write_bytes(b"x" * 1024)  # 1KB
        
        # Update and check
        size, delta, rate = monitor.update()
        assert size >= 1024, f"Size should be at least 1024 bytes, got {size}"
        print(f"  ✓ File size monitoring works: {size} bytes detected")
        
        # Update again (should show no change if file unchanged)
        time.sleep(0.6)
        size2, delta2, rate2 = monitor.update()
        assert delta2 == 0, "Delta should be 0 for unchanged file"
        print("  ✓ No-change detection works")
    
    print("✓ FileSizeMonitor tests passed\n")


def test_thread_progress_tracker():
    """Test ThreadProgressTracker functionality."""
    print("Testing ThreadProgressTracker...")
    
    with tempfile.TemporaryDirectory() as tmpdir:
        test_dir = Path(tmpdir) / "test_sample"
        test_dir.mkdir()
        
        tracker = ThreadProgressTracker(
            thread_id=1,
            run_id="SRR1234567",
            sample_dir=test_dir,
            update_interval=0.5,
            use_progress_bar=False,  # Disable bars for testing
        )
        
        assert not tracker.is_complete, "Tracker should not be complete initially"
        print("  ✓ Initial state check passed")
        
        # Create a test file
        test_file = test_dir / "test.fastq.gz"
        test_file.write_bytes(b"x" * 2048)  # 2KB
        
        # Update progress
        state = tracker.update()
        assert state["size_bytes"] >= 2048, "Should detect file size"
        assert state["elapsed_seconds"] >= 0, "Elapsed time should be non-negative"
        print(f"  ✓ Progress update works: {state['size_mb']:.3f} MB, {state['rate_mbps']:.2f} MB/s")
        
        # Mark complete
        tracker.mark_complete(success=True)
        assert tracker.is_complete, "Tracker should be marked complete"
        assert not tracker.failed, "Tracker should not be marked failed"
        print("  ✓ Completion marking works")
        
        tracker.close()
        print("  ✓ Cleanup works")
    
    print("✓ ThreadProgressTracker tests passed\n")


def test_download_progress_monitor():
    """Test DownloadProgressMonitor functionality."""
    print("Testing DownloadProgressMonitor...")
    
    with tempfile.TemporaryDirectory() as tmpdir:
        test_dir = Path(tmpdir)
        
        monitor = DownloadProgressMonitor(
            out_dir=test_dir,
            update_interval=0.5,
            use_progress_bars=False,  # Disable bars for testing
            show_summary=False,  # Disable summary for testing
        )
        
        # Register a thread
        monitor.register_thread(1, "SRR1234567")
        assert 1 in monitor.active_trackers, "Thread should be registered"
        print("  ✓ Thread registration works")
        
        # Create a test file
        sample_dir = test_dir / "getfastq" / "SRR1234567"
        sample_dir.mkdir(parents=True)
        test_file = sample_dir / "test.fastq.gz"
        test_file.write_bytes(b"x" * 4096)  # 4KB
        
        # Update thread
        state = monitor.update_thread(1)
        assert state is not None, "Should return state for registered thread"
        assert state["size_bytes"] >= 4096, "Should detect file size"
        print(f"  ✓ Thread update works: {state['size_mb']:.3f} MB")
        
        # Get summary
        summary = monitor.get_summary()
        assert summary["total_samples"] == 1, "Should have 1 sample"
        assert summary["active"] == 1, "Should have 1 active download"
        print(f"  ✓ Summary works: {summary['total_samples']} samples, {summary['active']} active")
        
        # Unregister thread
        monitor.unregister_thread(1, success=True)
        summary2 = monitor.get_summary()
        assert summary2["completed"] == 1, "Should have 1 completed"
        assert summary2["active"] == 0, "Should have 0 active"
        print("  ✓ Thread unregistration works")
        
        # Test monitoring
        monitor.start_monitoring()
        time.sleep(0.6)  # Wait for one update cycle
        monitor.stop_monitoring()
        print("  ✓ Monitoring start/stop works")
    
    print("✓ DownloadProgressMonitor tests passed\n")


def test_integration():
    """Test integration with multiple threads."""
    print("Testing multi-thread integration...")
    
    with tempfile.TemporaryDirectory() as tmpdir:
        test_dir = Path(tmpdir)
        
        monitor = DownloadProgressMonitor(
            out_dir=test_dir,
            update_interval=0.5,
            use_progress_bars=False,
            show_summary=False,
        )
        
        # Register multiple threads
        for i in range(3):
            run_id = f"SRR123456{i}"
            monitor.register_thread(i + 1, run_id)
            
            # Create test files
            sample_dir = test_dir / "getfastq" / run_id
            sample_dir.mkdir(parents=True)
            test_file = sample_dir / "test.fastq.gz"
            test_file.write_bytes(b"x" * (1024 * (i + 1)))  # Different sizes
        
        summary = monitor.get_summary()
        assert summary["total_samples"] == 3, "Should have 3 samples"
        assert summary["active"] == 3, "Should have 3 active downloads"
        print(f"  ✓ Multi-thread registration: {summary['total_samples']} samples")
        
        # Update all threads
        for i in range(3):
            state = monitor.update_thread(i + 1)
            assert state is not None, f"Thread {i+1} should have state"
        
        # Unregister all
        for i in range(3):
            monitor.unregister_thread(i + 1, success=True)
        
        final_summary = monitor.get_summary()
        assert final_summary["completed"] == 3, "Should have 3 completed"
        assert final_summary["active"] == 0, "Should have 0 active"
        print("  ✓ Multi-thread unregistration works")
        
        monitor.stop_monitoring()
    
    print("✓ Integration tests passed\n")


def main():
    """Run all tests."""
    print("=" * 80)
    print("COMPREHENSIVE PROGRESS TRACKING TESTS")
    print("=" * 80)
    print()
    
    try:
        test_file_size_monitor()
        test_thread_progress_tracker()
        test_download_progress_monitor()
        test_integration()
        
        print("=" * 80)
        print("✓ ALL TESTS PASSED")
        print("=" * 80)
        return 0
    except Exception as e:
        print(f"\n✗ TEST FAILED: {e}")
        import traceback
        traceback.print_exc()
        return 1


if __name__ == "__main__":
    sys.exit(main())

