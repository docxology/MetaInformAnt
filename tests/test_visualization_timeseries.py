"""Tests for time series visualization functions."""

from __future__ import annotations

from pathlib import Path
import numpy as np
import pandas as pd
import pytest

from metainformant.visualization.timeseries import (
    plot_time_series,
    plot_autocorrelation,
    plot_seasonal_decomposition,
    plot_forecast,
    plot_trend_analysis,
)

# Check for optional dependencies
try:
    from statsmodels.tsa.seasonal import seasonal_decompose
    HAS_STATSMODELS = True
except ImportError:
    HAS_STATSMODELS = False


class TestPlotTimeSeries:
    """Test plot_time_series function."""

    def test_basic_time_series_plot(self):
        """Test basic time series plot creation."""
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        # Create sample time series data
        dates = pd.date_range('2020-01-01', periods=100, freq='D')
        data = pd.DataFrame({
            'series1': np.random.randn(100).cumsum(),
            'series2': np.random.randn(100).cumsum() + 10
        }, index=dates)

        ax = plot_time_series(data)
        assert ax is not None
        assert len(ax.lines) == 2  # Two series plotted
        plt.close('all')

    def test_time_series_plot_single_series(self):
        """Test time series plot with single series."""
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        dates = pd.date_range('2020-01-01', periods=50, freq='D')
        data = pd.DataFrame({'value': np.sin(np.arange(50) * 0.1)}, index=dates)

        ax = plot_time_series(data)
        assert ax is not None
        assert len(ax.lines) == 1
        plt.close('all')

    def test_time_series_plot_with_output_path(self, tmp_path: Path):
        """Test time series plot with output path."""
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        dates = pd.date_range('2020-01-01', periods=30, freq='D')
        data = pd.DataFrame({'temp': np.random.randn(30)}, index=dates)
        output_path = tmp_path / "timeseries.png"

        ax = plot_time_series(data, output_path=output_path)
        assert ax is not None
        assert output_path.exists()
        plt.close('all')

    def test_time_series_plot_empty_data(self):
        """Test time series plot with empty data."""
        data = pd.DataFrame()

        with pytest.raises(ValueError, match="cannot be empty"):
            plot_time_series(data)


class TestPlotAutocorrelation:
    """Test plot_autocorrelation function."""

    def test_basic_autocorrelation_plot(self):
        """Test basic autocorrelation plot creation."""
        if not HAS_STATSMODELS:
            pytest.skip("statsmodels required for autocorrelation plotting")

        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        # Create sample time series
        dates = pd.date_range('2020-01-01', periods=100, freq='D')
        data = pd.Series(np.sin(np.arange(100) * 0.1) + np.random.randn(100) * 0.1,
                        index=dates)

        ax = plot_autocorrelation(data)
        assert ax is not None
        plt.close('all')

    def test_autocorrelation_plot_with_output_path(self, tmp_path: Path):
        """Test autocorrelation plot with output path."""
        if not HAS_STATSMODELS:
            pytest.skip("statsmodels required for autocorrelation plotting")

        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        data = pd.Series(np.random.randn(50))
        output_path = tmp_path / "autocorr.png"

        ax = plot_autocorrelation(data, output_path=output_path)
        assert ax is not None
        assert output_path.exists()
        plt.close('all')

    def test_autocorrelation_plot_empty_data(self):
        """Test autocorrelation plot with empty data."""
        if not HAS_STATSMODELS:
            pytest.skip("statsmodels required for autocorrelation plotting")

        data = pd.Series([], dtype=float)

        with pytest.raises(ValueError, match="cannot be empty"):
            plot_autocorrelation(data)

    def test_autocorrelation_plot_no_statsmodels(self):
        """Test autocorrelation plot when statsmodels is not available."""
        if HAS_STATSMODELS:
            pytest.skip("statsmodels is available")

        data = pd.Series(np.random.randn(50))

        with pytest.raises(ImportError, match="statsmodels required"):
            plot_autocorrelation(data)


class TestPlotSeasonalDecomposition:
    """Test plot_seasonal_decomposition function."""

    def test_basic_seasonal_decomposition_plot(self):
        """Test basic seasonal decomposition plot creation."""
        if not HAS_STATSMODELS:
            pytest.skip("statsmodels required for seasonal decomposition")

        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        # Create sample seasonal time series
        dates = pd.date_range('2020-01-01', periods=200, freq='D')
        seasonal = np.sin(2 * np.pi * np.arange(200) / 7)  # Weekly pattern
        trend = np.linspace(0, 10, 200)
        noise = np.random.randn(200) * 0.5
        data = pd.Series(seasonal + trend + noise, index=dates)

        ax = plot_seasonal_decomposition(data)
        assert ax is not None
        plt.close('all')

    def test_seasonal_decomposition_with_output_path(self, tmp_path: Path):
        """Test seasonal decomposition plot with output path."""
        if not HAS_STATSMODELS:
            pytest.skip("statsmodels required for seasonal decomposition")

        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        dates = pd.date_range('2020-01-01', periods=100, freq='D')
        data = pd.Series(np.sin(np.arange(100) * 0.1), index=dates)
        output_path = tmp_path / "decomp.png"

        ax = plot_seasonal_decomposition(data, output_path=output_path)
        assert ax is not None
        assert output_path.exists()
        plt.close('all')

    def test_seasonal_decomposition_multiplicative(self):
        """Test seasonal decomposition with multiplicative model."""
        if not HAS_STATSMODELS:
            pytest.skip("statsmodels required for seasonal decomposition")

        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        dates = pd.date_range('2020-01-01', periods=100, freq='D')
        data = pd.Series(np.abs(np.sin(np.arange(100) * 0.1)) + 1, index=dates)

        ax = plot_seasonal_decomposition(data, model='multiplicative')
        assert ax is not None
        plt.close('all')

    def test_seasonal_decomposition_empty_data(self):
        """Test seasonal decomposition plot with empty data."""
        if not HAS_STATSMODELS:
            pytest.skip("statsmodels required for seasonal decomposition")

        data = pd.Series([], dtype=float)

        with pytest.raises(ValueError, match="cannot be empty"):
            plot_seasonal_decomposition(data)

    def test_seasonal_decomposition_no_statsmodels(self):
        """Test seasonal decomposition when statsmodels is not available."""
        if HAS_STATSMODELS:
            pytest.skip("statsmodels is available")

        data = pd.Series(np.random.randn(50))

        with pytest.raises(ImportError, match="statsmodels required"):
            plot_seasonal_decomposition(data)


class TestPlotForecast:
    """Test plot_forecast function."""

    def test_basic_forecast_plot(self):
        """Test basic forecast plot creation."""
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        # Create historical and forecast data
        dates_hist = pd.date_range('2020-01-01', periods=50, freq='D')
        dates_forecast = pd.date_range('2020-02-20', periods=30, freq='D')

        historical = pd.Series(np.sin(np.arange(50) * 0.1), index=dates_hist)
        forecast = pd.Series(np.sin(np.arange(50, 80) * 0.1), index=dates_forecast)

        ax = plot_forecast(historical, forecast)
        assert ax is not None
        assert len(ax.lines) == 2  # Historical and forecast lines
        plt.close('all')

    def test_forecast_plot_with_confidence_intervals(self):
        """Test forecast plot with confidence intervals."""
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        dates_hist = pd.date_range('2020-01-01', periods=30, freq='D')
        dates_forecast = pd.date_range('2020-01-31', periods=10, freq='D')

        historical = pd.Series(np.random.randn(30), index=dates_hist)
        forecast = pd.Series(np.random.randn(10), index=dates_forecast)

        # Create confidence intervals
        ci = pd.DataFrame({
            'lower': forecast - 0.5,
            'upper': forecast + 0.5
        }, index=dates_forecast)

        ax = plot_forecast(historical, forecast, confidence_intervals=ci)
        assert ax is not None
        plt.close('all')

    def test_forecast_plot_with_output_path(self, tmp_path: Path):
        """Test forecast plot with output path."""
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        dates_hist = pd.date_range('2020-01-01', periods=20, freq='D')
        dates_forecast = pd.date_range('2020-01-21', periods=10, freq='D')

        historical = pd.Series(np.arange(20), index=dates_hist)
        forecast = pd.Series(np.arange(20, 30), index=dates_forecast)
        output_path = tmp_path / "forecast.png"

        ax = plot_forecast(historical, forecast, output_path=output_path)
        assert ax is not None
        assert output_path.exists()
        plt.close('all')

    def test_forecast_plot_empty_data(self):
        """Test forecast plot with empty data."""
        historical = pd.Series([], dtype=float)
        forecast = pd.Series([1, 2, 3])

        with pytest.raises(ValueError, match="cannot be empty"):
            plot_forecast(historical, forecast)

    def test_forecast_plot_empty_forecast(self):
        """Test forecast plot with empty forecast."""
        historical = pd.Series([1, 2, 3])
        forecast = pd.Series([], dtype=float)

        with pytest.raises(ValueError, match="cannot be empty"):
            plot_forecast(historical, forecast)


class TestPlotTrendAnalysis:
    """Test plot_trend_analysis function."""

    def test_basic_trend_analysis_plot(self):
        """Test basic trend analysis plot creation."""
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        # Create time series with clear trend
        dates = pd.date_range('2020-01-01', periods=100, freq='D')
        trend = np.linspace(0, 10, 100)
        noise = np.random.randn(100) * 0.5
        data = pd.Series(trend + noise, index=dates)

        ax = plot_trend_analysis(data)
        assert ax is not None
        assert len(ax.lines) == 2  # Original data and trend line
        plt.close('all')

    def test_trend_analysis_with_output_path(self, tmp_path: Path):
        """Test trend analysis plot with output path."""
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        dates = pd.date_range('2020-01-01', periods=50, freq='D')
        data = pd.Series(np.random.randn(50).cumsum(), index=dates)
        output_path = tmp_path / "trend.png"

        ax = plot_trend_analysis(data, output_path=output_path)
        assert ax is not None
        assert output_path.exists()
        plt.close('all')

    def test_trend_analysis_empty_data(self):
        """Test trend analysis plot with empty data."""
        data = pd.Series([], dtype=float)

        with pytest.raises(ValueError, match="cannot be empty"):
            plot_trend_analysis(data)






