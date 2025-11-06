#!/usr/bin/env python3
"""Comprehensive test suite for download progress tracking with realistic scenarios."""

import sys
import tempfile
import threading
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from metainformant.rna.steps.download_progress import (
    DownloadProgressMonitor,
    ThreadProgressTracker,
    FileSizeMonitor,
)


def test_concurrent_file_updates():
    """Test that file size monitoring works with concurrent file updates."""
    print("Testing concurrent file updates...")
    
    with tempfile.TemporaryDirectory() as tmpdir:
        test_dir = Path(tmpdir) / "test_sample"
        test_dir.mkdir()
        
        monitor = FileSizeMonitor(test_dir, update_interval=0.1)
        
        def write_file(thread_id: int, iterations: int):
            """Simulate file writing."""
            for i in range(iterations):
                test_file = test_dir / f"file_{thread_id}_{i}.fastq.gz"
                test_file.write_bytes(b"x" * 1024)  # 1KB per file
                time.sleep(0.05)
        
        # Start multiple threads writing files
        threads = []
        for i in range(3):
            t = threading.Thread(target=write_file, args=(i, 5))
            t.start()
            threads.append(t)
        
        # Monitor while files are being written
        updates = []
        for _ in range(10):
            size, delta, rate = monitor.update()
            updates.append((size, delta, rate))
            time.sleep(0.15)
        
        # Wait for all threads
        for t in threads:
            t.join()
        
        # Final update
        final_size, _, _ = monitor.update()
        
        assert final_size >= 15 * 1024, f"Should have at least 15KB (15 files), got {final_size}"
        assert any(delta > 0 for _, delta, _ in updates), "Should detect file size changes"
        print(f"  ✓ Detected concurrent file updates: {final_size} bytes total")
        print(f"  ✓ Rate calculation works: {updates[-1][2]:.2f} MB/s")
    
    print("✓ Concurrent file updates test passed\n")


def test_thread_safety():
    """Test thread safety of DownloadProgressMonitor."""
    print("Testing thread safety...")
    
    with tempfile.TemporaryDirectory() as tmpdir:
        test_dir = Path(tmpdir)
        
        monitor = DownloadProgressMonitor(
            out_dir=test_dir,
            update_interval=0.1,
            use_progress_bars=False,
            show_summary=False,
        )
        
        # Register multiple threads concurrently
        def register_and_update(thread_id: int, run_id: str):
            """Register thread and update progress."""
            monitor.register_thread(thread_id, run_id)
            
            # Create sample directory and file
            sample_dir = test_dir / "getfastq" / run_id
            sample_dir.mkdir(parents=True)
            
            # Simulate file growth
            for i in range(5):
                test_file = sample_dir / f"chunk_{i}.fastq.gz"
                test_file.write_bytes(b"x" * (1024 * (i + 1)))
                monitor.update_thread(thread_id)
                time.sleep(0.05)
            
            monitor.unregister_thread(thread_id, success=True)
        
        # Start multiple threads
        threads = []
        for i in range(5):
            t = threading.Thread(target=register_and_update, args=(i + 1, f"SRR123456{i}"))
            t.start()
            threads.append(t)
        
        # Wait for all threads
        for t in threads:
            t.join()
        
        summary = monitor.get_summary()
        assert summary["total_samples"] == 5, "Should have 5 samples"
        assert summary["completed"] == 5, "Should have 5 completed"
        assert summary["active"] == 0, "Should have 0 active"
        print(f"  ✓ Thread-safe operations: {summary['total_samples']} samples processed")
        
        monitor.stop_monitoring()
    
    print("✓ Thread safety test passed\n")


def test_error_handling():
    """Test error handling in progress monitoring."""
    print("Testing error handling...")
    
    with tempfile.TemporaryDirectory() as tmpdir:
        test_dir = Path(tmpdir)
        
        # Test with non-existent directory
        monitor = FileSizeMonitor(Path("/nonexistent/directory"), update_interval=0.1)
        size = monitor.get_size()
        assert size == 0, "Should return 0 for non-existent directory"
        print("  ✓ Handles non-existent directories gracefully")
        
        # Test with deleted directory
        temp_sample_dir = test_dir / "temp_sample"
        temp_sample_dir.mkdir()
        monitor2 = FileSizeMonitor(temp_sample_dir, update_interval=0.1)
        
        # Create and delete file
        test_file = temp_sample_dir / "test.fastq.gz"
        test_file.write_bytes(b"x" * 1024)
        monitor2.update()
        
        test_file.unlink()
        size2, _, _ = monitor2.update()
        assert size2 == 0, "Should handle deleted files"
        print("  ✓ Handles deleted files gracefully")
        
        # Test monitor with invalid thread IDs
        monitor3 = DownloadProgressMonitor(
            out_dir=test_dir,
            update_interval=0.1,
            use_progress_bars=False,
            show_summary=False,
        )
        
        # Try to update non-existent thread
        state = monitor3.update_thread(999)
        assert state is None, "Should return None for non-existent thread"
        print("  ✓ Handles invalid thread IDs gracefully")
        
        # Try to unregister non-existent thread (should not crash)
        try:
            monitor3.unregister_thread(999, success=True)
            print("  ✓ Handles unregistering non-existent threads gracefully")
        except Exception as e:
            print(f"  ✗ Error unregistering non-existent thread: {e}")
            raise
        
        monitor3.stop_monitoring()
    
    print("✓ Error handling test passed\n")


def test_realistic_download_simulation():
    """Simulate a realistic download scenario with multiple samples."""
    print("Testing realistic download simulation...")
    
    with tempfile.TemporaryDirectory() as tmpdir:
        test_dir = Path(tmpdir)
        
        monitor = DownloadProgressMonitor(
            out_dir=test_dir,
            update_interval=0.2,
            use_progress_bars=False,
            show_summary=True,
        )
        
        monitor.start_monitoring()
        
        # Simulate 3 samples downloading at different rates
        samples = [
            ("SRR1234567", 0.1, 10),  # Fast download: 0.1s per chunk, 10 chunks
            ("SRR1234568", 0.2, 8),   # Medium download: 0.2s per chunk, 8 chunks
            ("SRR1234569", 0.3, 6),   # Slow download: 0.3s per chunk, 6 chunks
        ]
        
        def simulate_download(run_id: str, delay: float, chunks: int, thread_id: int):
            """Simulate a download."""
            monitor.register_thread(thread_id, run_id)
            
            sample_dir = test_dir / "getfastq" / run_id
            sample_dir.mkdir(parents=True)
            
            for i in range(chunks):
                # Simulate chunked download
                chunk_file = sample_dir / f"chunk_{i}.fastq.gz"
                chunk_file.write_bytes(b"x" * (1024 * 100))  # 100KB per chunk
                time.sleep(delay)
            
            # Final file
            final_file = sample_dir / f"{run_id}_1.fastq.gz"
            final_file.write_bytes(b"x" * (1024 * 50))  # 50KB final file
            
            monitor.unregister_thread(thread_id, success=True)
        
        # Start downloads
        threads = []
        for idx, (run_id, delay, chunks) in enumerate(samples, 1):
            t = threading.Thread(target=simulate_download, args=(run_id, delay, chunks, idx))
            t.start()
            threads.append(t)
        
        # Wait for all downloads
        for t in threads:
            t.join()
        
        # Give monitor time to update
        time.sleep(0.5)
        
        summary = monitor.get_summary()
        assert summary["total_samples"] == 3, "Should have 3 samples"
        assert summary["completed"] == 3, "Should have 3 completed"
        assert summary["total_size_mb"] > 0, "Should have downloaded data"
        print(f"  ✓ Simulated download: {summary['total_samples']} samples")
        print(f"  ✓ Total size: {summary['total_size_mb']:.2f} MB")
        print(f"  ✓ Average rate: {summary['avg_rate_mbps']:.2f} MB/s")
        
        monitor.stop_monitoring()
    
    print("✓ Realistic download simulation test passed\n")


def test_progress_bar_fallback():
    """Test that progress bars gracefully fall back to text if tqdm unavailable."""
    print("Testing progress bar fallback...")
    
    with tempfile.TemporaryDirectory() as tmpdir:
        test_dir = Path(tmpdir) / "test_sample"
        test_dir.mkdir()
        
        # Create tracker with progress bar disabled (simulating tqdm unavailable)
        tracker = ThreadProgressTracker(
            thread_id=1,
            run_id="SRR1234567",
            sample_dir=test_dir,
            update_interval=0.1,
            use_progress_bar=False,  # Simulate no tqdm
        )
        
        # Should still work without progress bars
        test_file = test_dir / "test.fastq.gz"
        test_file.write_bytes(b"x" * 2048)
        
        state = tracker.update()
        assert state["size_bytes"] >= 2048, "Should work without progress bars"
        print("  ✓ Works without progress bars (text mode)")
        
        tracker.mark_complete(success=True)
        tracker.close()
    
    print("✓ Progress bar fallback test passed\n")


def test_configuration_variations():
    """Test different configuration combinations."""
    print("Testing configuration variations...")
    
    with tempfile.TemporaryDirectory() as tmpdir:
        test_dir = Path(tmpdir)
        
        configs = [
            {"use_progress_bars": True, "show_summary": False},
            {"use_progress_bars": False, "show_summary": True},
            {"use_progress_bars": False, "show_summary": False},
            {"update_interval": 1.0},
            {"update_interval": 5.0},
        ]
        
        for i, config in enumerate(configs):
            monitor = DownloadProgressMonitor(
                out_dir=test_dir,
                update_interval=config.get("update_interval", 2.0),
                use_progress_bars=config.get("use_progress_bars", True),
                show_summary=config.get("show_summary", False),
            )
            
            assert monitor.update_interval == config.get("update_interval", 2.0)
            assert monitor.use_progress_bars == config.get("use_progress_bars", True)
            assert monitor.show_summary == config.get("show_summary", False)
            
            monitor.stop_monitoring()
        
        print(f"  ✓ Tested {len(configs)} different configurations")
    
    print("✓ Configuration variations test passed\n")


def test_large_file_handling():
    """Test handling of large files and many files."""
    print("Testing large file handling...")
    
    with tempfile.TemporaryDirectory() as tmpdir:
        test_dir = Path(tmpdir) / "test_sample"
        test_dir.mkdir()
        
        monitor = FileSizeMonitor(test_dir, update_interval=0.1)
        
        # Create many small files
        for i in range(100):
            test_file = test_dir / f"file_{i}.fastq.gz"
            test_file.write_bytes(b"x" * 1024)  # 1KB each
        
        size1 = monitor.get_size()
        assert size1 >= 100 * 1024, f"Should detect 100KB, got {size1}"
        print(f"  ✓ Handles many files: {size1} bytes from 100 files")
        
        # Create a large file
        large_file = test_dir / "large.fastq.gz"
        large_file.write_bytes(b"x" * (10 * 1024 * 1024))  # 10MB
        
        size2 = monitor.get_size()
        assert size2 >= 10 * 1024 * 1024, f"Should detect 10MB, got {size2}"
        print(f"  ✓ Handles large files: {size2 / (1024*1024):.1f} MB")
    
    print("✓ Large file handling test passed\n")


def main():
    """Run all comprehensive tests."""
    print("=" * 80)
    print("COMPREHENSIVE PROGRESS TRACKING TEST SUITE")
    print("=" * 80)
    print()
    
    tests = [
        test_concurrent_file_updates,
        test_thread_safety,
        test_error_handling,
        test_realistic_download_simulation,
        test_progress_bar_fallback,
        test_configuration_variations,
        test_large_file_handling,
    ]
    
    passed = 0
    failed = 0
    
    for test_func in tests:
        try:
            test_func()
            passed += 1
        except Exception as e:
            print(f"\n✗ {test_func.__name__} FAILED: {e}")
            import traceback
            traceback.print_exc()
            failed += 1
            print()
    
    print("=" * 80)
    print(f"TEST RESULTS: {passed} passed, {failed} failed")
    print("=" * 80)
    
    if failed == 0:
        print("✓ ALL COMPREHENSIVE TESTS PASSED")
        return 0
    else:
        print(f"✗ {failed} TEST(S) FAILED")
        return 1


if __name__ == "__main__":
    sys.exit(main())

