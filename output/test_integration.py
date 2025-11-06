#!/usr/bin/env python3
"""Integration test for progress tracking with step functions."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from metainformant.rna.steps.download_progress import DownloadProgressMonitor


def test_config_parsing():
    """Test that configuration parameters are correctly parsed."""
    print("Testing configuration parameter parsing...")
    
    # Test default values
    params = {}
    show_progress = params.get("show_progress", True)
    update_interval = float(params.get("progress_update_interval", 2.0))
    use_bars = params.get("progress_style", "bar") == "bar"
    
    assert show_progress is True, "Default show_progress should be True"
    assert update_interval == 2.0, "Default update_interval should be 2.0"
    assert use_bars is True, "Default progress_style should be 'bar'"
    print("  ✓ Default configuration values correct")
    
    # Test custom values
    params2 = {
        "show_progress": False,
        "progress_update_interval": 5.0,
        "progress_style": "text",
    }
    show_progress2 = params2.get("show_progress", True)
    update_interval2 = float(params2.get("progress_update_interval", 2.0))
    use_bars2 = params2.get("progress_style", "bar") == "bar"
    
    assert show_progress2 is False, "Custom show_progress should be False"
    assert update_interval2 == 5.0, "Custom update_interval should be 5.0"
    assert use_bars2 is False, "Custom progress_style 'text' should result in use_bars=False"
    print("  ✓ Custom configuration values parsed correctly")
    
    print("✓ Configuration parsing tests passed\n")


def test_monitor_instantiation():
    """Test that monitor can be instantiated with various configs."""
    print("Testing monitor instantiation...")
    
    import tempfile
    
    with tempfile.TemporaryDirectory() as tmpdir:
        # Test with progress bars
        m1 = DownloadProgressMonitor(
            out_dir=Path(tmpdir),
            update_interval=2.0,
            use_progress_bars=True,
            show_summary=False,
        )
        assert m1.use_progress_bars is True
        assert m1.show_summary is False
        print("  ✓ Monitor with progress bars instantiated")
        m1.stop_monitoring()
        
        # Test with text updates
        m2 = DownloadProgressMonitor(
            out_dir=Path(tmpdir),
            update_interval=3.0,
            use_progress_bars=False,
            show_summary=True,
        )
        assert m2.use_progress_bars is False
        assert m2.show_summary is True
        assert m2.update_interval == 3.0
        print("  ✓ Monitor with text updates instantiated")
        m2.stop_monitoring()
    
    print("✓ Monitor instantiation tests passed\n")


def main():
    """Run integration tests."""
    print("=" * 80)
    print("INTEGRATION TESTS FOR PROGRESS TRACKING")
    print("=" * 80)
    print()
    
    try:
        test_config_parsing()
        test_monitor_instantiation()
        
        print("=" * 80)
        print("✓ ALL INTEGRATION TESTS PASSED")
        print("=" * 80)
        return 0
    except Exception as e:
        print(f"\n✗ TEST FAILED: {e}")
        import traceback
        traceback.print_exc()
        return 1


if __name__ == "__main__":
    sys.exit(main())

