from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from metainformant.core import io as core_io
from .model import GenerationResult, GenerationsResult, colors, lin_phi_bar, lin_phi_inv


def _finalize_axes(ax: plt.Axes, *, show_ticks: bool = False) -> None:
    """Apply consistent styling for accessibility and readability."""
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_linewidth(2)
    ax.spines['bottom'].set_color('#333')
    ax.spines['left'].set_linewidth(2)
    ax.spines['left'].set_color('#333')
    ax.grid(False)

    # Set larger, accessible font sizes
    ax.tick_params(axis='both', which='major', labelsize=14, width=2, length=6)
    ax.xaxis.label.set_size(16)
    ax.yaxis.label.set_size(16)

    if not show_ticks:
        ax.set_xticks([])
        ax.set_yticks([])


def display_s_vs_q(r: GenerationResult, *, sample_size: int = 600, dest: str | Path | None = None) -> Path | None:
    """Create accessible visualization of structural trait vs quality relationship."""
    indices = np.random.choice(r.s.shape[0], sample_size, replace=False)
    s_sample = r.s[indices]
    q_sample = r.q[indices]

    fig, ax = plt.subplots(figsize=(12, 10), dpi=200, facecolor='white')
    _finalize_axes(ax, show_ticks=True)

    # Add subtle grid for better readability
    ax.grid(True, alpha=0.3, linestyle='--', linewidth=0.5)

    # Plot the data points with better contrast and size
    plt.scatter(s_sample, q_sample, alpha=0.6, s=20, color=colors['green'],
                edgecolor='darkgreen', linewidth=0.5, label=r"$q$ (quality signals)")

    # Plot the true relationship line
    s_line = np.linspace(np.min(s_sample), np.max(s_sample), 100)
    plt.plot(s_line, lin_phi_bar(s_line), color=colors['blue'], linewidth=3,
             linestyle='-', label=r"$\bar{\phi}(s)$ (true relationship)")

    # Improve labels with larger fonts
    plt.xlabel(r"Structural Trait $s$", fontsize=18, fontweight='bold')
    plt.ylabel(r"Quality Signal $q$", fontsize=18, fontweight='bold')

    # Create informative title and legend
    plt.title(f"Structural Trait vs Quality Signal\nSignal Fidelity: $\\hat{{s}}$ = {r.s_hat}",
              fontsize=20, fontweight='bold', pad=30, loc='center')

    # Position legend in better location with larger font
    legend = plt.legend(fontsize=16, loc='upper left', framealpha=0.9,
                       bbox_to_anchor=(0.02, 0.98))
    legend.get_frame().set_linewidth(2)

    # Add correlation information as text annotation
    correlation_text = f'Correlation: ρ(s,q) = {r.rho_sq:.3f}'
    ax.text(0.02, 0.02, correlation_text, transform=ax.transAxes,
            fontsize=14, verticalalignment='bottom',
            bbox=dict(boxstyle='round,pad=0.3', facecolor='lightblue', alpha=0.8))

    out_path: Path | None = None
    if dest is not None:
        out_path = Path(dest)
        core_io.ensure_directory(out_path.parent)
        plt.savefig(out_path, format='png', dpi=200, bbox_inches='tight',
                   facecolor='white', edgecolor='none')
    plt.close(fig)
    return out_path


def display_sq_vs_w(r: GenerationResult, *, sample_size: int = 600, dest: str | Path | None = None) -> Path | None:
    """Create accessible visualization of traits vs fitness relationships."""
    indices = np.random.choice(r.s.shape[0], sample_size, replace=False)
    s_sample = r.s[indices]
    q_sample = r.q[indices]
    w_sample = r.w[indices]

    fig, ax = plt.subplots(figsize=(12, 10), dpi=200, facecolor='white')
    _finalize_axes(ax, show_ticks=True)

    # Add subtle grid
    ax.grid(True, alpha=0.3, linestyle='--', linewidth=0.5)

    s_from_q = lin_phi_inv(q_sample)
    coeffs_s = np.polyfit(s_sample, w_sample, 1)
    fit_line_s = np.poly1d(coeffs_s)
    coeffs_q = np.polyfit(s_from_q, w_sample, 1)
    fit_line_q = np.poly1d(coeffs_q)

    x_vals = np.linspace(
        min(np.min(s_sample), np.min(s_from_q)),
        max(np.max(s_sample), np.max(s_from_q)),
        100,
    )

    # Plot data points with better visual distinction
    plt.scatter(s_sample, w_sample, alpha=0.7, s=25, color=colors['blue'],
                edgecolor='navy', linewidth=0.5, label=r"$s$ (structural traits)")
    plt.scatter(s_from_q, w_sample, alpha=0.7, s=25, color=colors['red'],
                edgecolor='darkred', linewidth=0.5, label=r"$q$ (quality signals)")

    # Plot regression lines with better styling
    plt.plot(x_vals, fit_line_s(x_vals), color=colors['green'], linewidth=3,
             linestyle='-', alpha=0.8, label=".2f")
    plt.plot(x_vals, fit_line_q(x_vals), color=colors['yellow'], linewidth=3,
             linestyle='--', alpha=0.8, label=".2f")

    # Improved labels and title
    plt.xlabel(r"Trait Value", fontsize=18, fontweight='bold')
    plt.ylabel(r"Fitness $w$", fontsize=18, fontweight='bold')
    plt.title(f"Selection on Correlated Traits\nSignal Fidelity: $\\hat{{s}}$ = {r.s_hat}",
              fontsize=20, fontweight='bold', pad=30)

    # Enhanced legend
    legend = plt.legend(fontsize=16, loc='upper right', framealpha=0.9,
                       bbox_to_anchor=(0.98, 0.98))
    legend.get_frame().set_linewidth(2)

    # Add selection differential information
    selection_info = (f'Selection on Structure: {r.rho_ws:.3f}\n'
                     f'Selection on Signal: {r.rho_wq:.3f}\n'
                     f'Relative Selection: {r.ratio_q:.3f}')
    ax.text(0.02, 0.02, selection_info, transform=ax.transAxes,
            fontsize=14, verticalalignment='bottom',
            bbox=dict(boxstyle='round,pad=0.5', facecolor='lightyellow', alpha=0.9))

    out_path: Path | None = None
    if dest is not None:
        out_path = Path(dest)
        core_io.ensure_directory(out_path.parent)
        plt.savefig(out_path, format='png', dpi=200, bbox_inches='tight',
                   facecolor='white', edgecolor='none')
    plt.close(fig)
    return out_path


def display_gr(gr: GenerationsResult, *, dest: str | Path | None = None) -> Path | None:
    """Create accessible visualization of trait evolution over generations."""
    start = 0
    step = (gr.generations // 20) + 1

    # Prepare data for plotting with better normalization
    s_means_p = gr.s_means[start::step].copy()
    s_stds_p = gr.s_stds[start::step].copy()

    q_means_p = lin_phi_inv(gr.q_means[start::step]).copy()
    q_stds_p = gr.q_stds[start::step].copy()

    w_means_p = gr.w_means[start::step].copy()
    w_stds_p = gr.w_stds[start::step].copy()

    # Normalize to 0-1 range for better comparison
    s_min, s_max = np.min(s_means_p), np.max(s_means_p)
    q_min, q_max = np.min(q_means_p), np.max(q_means_p)
    w_min, w_max = np.min(w_means_p), np.max(w_means_p)

    s_means_norm = (s_means_p - s_min) / (s_max - s_min) if s_max > s_min else s_means_p
    q_means_norm = (q_means_p - q_min) / (q_max - q_min) if q_max > q_min else q_means_p
    w_means_norm = (w_means_p - w_min) / (w_max - w_min) if w_max > w_min else w_means_p

    s_stds_norm = s_stds_p / (s_max - s_min) if s_max > s_min else s_stds_p
    q_stds_norm = q_stds_p / (q_max - q_min) if q_max > q_min else q_stds_p
    w_stds_norm = w_stds_p / (w_max - w_min) if w_max > w_min else w_stds_p

    fig, ax = plt.subplots(figsize=(14, 10), dpi=200, facecolor='white')
    _finalize_axes(ax, show_ticks=True)

    # Add grid for better readability
    ax.grid(True, alpha=0.3, linestyle='--', linewidth=0.5)

    # Create generation numbers for x-axis
    generations = [i * step for i in range(len(s_means_norm))]

    # Plot trait evolution lines with better styling
    plt.plot(generations, s_means_norm, color=colors['blue'], linewidth=3,
             marker='o', markersize=6, label=r"$s$ (structural traits)")
    plt.fill_between(generations,
                     s_means_norm - s_stds_norm,
                     s_means_norm + s_stds_norm,
                     color=colors['blue'], alpha=0.2)

    plt.plot(generations, q_means_norm, color=colors['red'], linewidth=3,
             marker='s', markersize=6, label=r"$q$ (quality signals)")
    plt.fill_between(generations,
                     q_means_norm - q_stds_norm,
                     q_means_norm + q_stds_norm,
                     color=colors['red'], alpha=0.2)

    plt.plot(generations, w_means_norm, color=colors['yellow'], linewidth=3,
             marker='^', markersize=6, label=r"$w$ (fitness)")
    plt.fill_between(generations,
                     w_means_norm - w_stds_norm,
                     w_means_norm + w_stds_norm,
                     color=colors['yellow'], alpha=0.2)

    # Enhanced labels and title
    plt.xlabel(r"Generation", fontsize=18, fontweight='bold')
    plt.ylabel(r"Normalized Trait Values", fontsize=18, fontweight='bold')
    plt.title(f"Trait Evolution Over {gr.generations} Generations\n"
              f"Signal Fidelity: $\\hat{{s}}$ = {gr.last_result.s_hat}",
              fontsize=20, fontweight='bold', pad=30)

    # Enhanced legend
    legend = plt.legend(fontsize=16, loc='center right', framealpha=0.9,
                       bbox_to_anchor=(1.02, 0.5))
    legend.get_frame().set_linewidth(2)

    # Add final statistics as text annotation
    final_stats = (f'Final Statistics:\n'
                  f'ρ(s,q) = {gr.rho_sq_means:.3f}\n'
                  f'QSC = {gr.last_result.QSC:.3f}\n'
                  f'QSC ratio = {gr.last_result.QSC_ratio:.3f}')
    ax.text(0.02, 0.98, final_stats, transform=ax.transAxes,
            fontsize=14, verticalalignment='top',
            bbox=dict(boxstyle='round,pad=0.5', facecolor='lightgray', alpha=0.9))

    out_path: Path | None = None
    if dest is not None:
        out_path = Path(dest)
        core_io.ensure_directory(out_path.parent)
        plt.savefig(out_path, format='png', dpi=200, bbox_inches='tight',
                   facecolor='white', edgecolor='none')
    plt.close(fig)
    return out_path


def display_generation(gr: GenerationsResult, *, dest: str | Path | None = None) -> Path | None:
    """Display generation statistics and evolution plot."""
    print(f"rho_sq_mean: {gr.rho_sq_means:.5f}")
    print("")
    print(f"q rho   : {gr.last_result.rho_wq:.5f}")
    print(f"q mean  : {np.mean(gr.last_result.q):.5f}")
    print(f"q CV    : {np.std(gr.last_result.q) / np.mean(gr.last_result.q):.5f}")
    print(f"q ratio : {gr.last_result.ratio_q:.5f}")
    print("")
    print(f"s rho   : {gr.last_result.rho_ws:.5f}")
    print(f"s mean  : {np.mean(gr.last_result.s):.5f}")
    print(f"s CV    : {np.std(gr.last_result.s) / np.mean(gr.last_result.s):.5f}")
    print(f"s ratio : {gr.last_result.ratio_s:.5f}")
    print("")
    print(f"QSC      : {gr.last_result.QSC:.5f}")
    print(f"QSC ratio: {gr.last_result.QSC_ratio:.5f}")
    return display_gr(gr, dest=dest)


def create_composite_abstract(experiment_results: dict, *, dest: str | Path | None = None) -> Path | None:
    """Create a composite graphical abstract combining all experiment results.

    Args:
        experiment_results: Dictionary with keys as experiment names and values as result objects
        dest: Destination path for the composite figure

    Returns:
        Path to saved composite figure or None if dest is None
    """
    # Create a large figure with subplots
    fig = plt.figure(figsize=(20, 16), dpi=300, facecolor='white')

    # Main title
    fig.suptitle("Natural Selection Experiments:\nSignal Processing in Evolutionary Biology",
                fontsize=28, fontweight='bold', y=0.98)

    # Create grid layout
    gs = fig.add_gridspec(3, 4, hspace=0.3, wspace=0.3,
                          top=0.9, bottom=0.1, left=0.1, right=0.95)

    # Panel A: Signal mapping relationship
    ax1 = fig.add_subplot(gs[0, 0:2])
    if 'signal_mapping' in experiment_results:
        r = experiment_results['signal_mapping']
        indices = np.random.choice(r.s.shape[0], 300, replace=False)
        s_sample = r.s[indices]
        q_sample = r.q[indices]

        ax1.scatter(s_sample, q_sample, alpha=0.6, s=30, color=colors['green'],
                   edgecolor='darkgreen', linewidth=0.5)
        s_line = np.linspace(np.min(s_sample), np.max(s_sample), 100)
        ax1.plot(s_line, lin_phi_bar(s_line), color=colors['blue'], linewidth=3)

        ax1.set_xlabel('Structural Trait (s)', fontsize=16, fontweight='bold')
        ax1.set_ylabel('Quality Signal (q)', fontsize=16, fontweight='bold')
        ax1.set_title('A. Signal Mapping Function', fontsize=18, fontweight='bold', pad=20)
        ax1.grid(True, alpha=0.3)
        ax1.text(0.02, 0.98, f'ρ(s,q) = {r.rho_sq:.3f}\nSignal fidelity = {r.s_hat}',
                transform=ax1.transAxes, fontsize=14, verticalalignment='top',
                bbox=dict(boxstyle='round,pad=0.3', facecolor='lightblue', alpha=0.8))

    # Panel B: Selection on correlated traits
    ax2 = fig.add_subplot(gs[0, 2:4])
    if 'selection_correlation' in experiment_results:
        r = experiment_results['selection_correlation']
        indices = np.random.choice(r.s.shape[0], 300, replace=False)
        s_sample = r.s[indices]
        q_sample = r.q[indices]
        w_sample = r.w[indices]

        s_from_q = lin_phi_inv(q_sample)
        coeffs_s = np.polyfit(s_sample, w_sample, 1)
        coeffs_q = np.polyfit(s_from_q, w_sample, 1)
        x_vals = np.linspace(min(np.min(s_sample), np.min(s_from_q)),
                           max(np.max(s_sample), np.max(s_from_q)), 100)

        ax2.scatter(s_sample, w_sample, alpha=0.7, s=30, color=colors['blue'],
                   edgecolor='navy', linewidth=0.5, label='Structural (s)')
        ax2.scatter(s_from_q, w_sample, alpha=0.7, s=30, color=colors['red'],
                   edgecolor='darkred', linewidth=0.5, label='Quality (q)')
        ax2.plot(x_vals, np.poly1d(coeffs_s)(x_vals), color=colors['green'], linewidth=3, alpha=0.8)
        ax2.plot(x_vals, np.poly1d(coeffs_q)(x_vals), color=colors['yellow'], linewidth=3, alpha=0.8)

        ax2.set_xlabel('Trait Value', fontsize=16, fontweight='bold')
        ax2.set_ylabel('Fitness (w)', fontsize=16, fontweight='bold')
        ax2.set_title('B. Selection on Correlated Traits', fontsize=18, fontweight='bold', pad=20)
        ax2.legend(fontsize=14)
        ax2.grid(True, alpha=0.3)

    # Panel C: Trait evolution over generations (rebound)
    ax3 = fig.add_subplot(gs[1, 0:2])
    if 'rebound_evolution' in experiment_results:
        gr = experiment_results['rebound_evolution']
        start, step = 0, (gr.generations // 15) + 1
        generations = [i * step for i in range(len(gr.s_means[start::step]))]

        # Normalize data
        s_means = gr.s_means[start::step]
        q_means = lin_phi_inv(gr.q_means[start::step])
        w_means = gr.w_means[start::step]

        s_norm = (s_means - np.min(s_means)) / (np.max(s_means) - np.min(s_means))
        q_norm = (q_means - np.min(q_means)) / (np.max(q_means) - np.min(q_means))
        w_norm = (w_means - np.min(w_means)) / (np.max(w_means) - np.min(w_means))

        ax3.plot(generations, s_norm, color=colors['blue'], linewidth=3, marker='o',
                markersize=4, label='Structural (s)')
        ax3.plot(generations, q_norm, color=colors['red'], linewidth=3, marker='s',
                markersize=4, label='Quality (q)')
        ax3.plot(generations, w_norm, color=colors['yellow'], linewidth=3, marker='^',
                markersize=4, label='Fitness (w)')

        ax3.set_xlabel('Generation', fontsize=16, fontweight='bold')
        ax3.set_ylabel('Normalized Value', fontsize=16, fontweight='bold')
        ax3.set_title('C. Rebound Selection Dynamics', fontsize=18, fontweight='bold', pad=20)
        ax3.legend(fontsize=14)
        ax3.grid(True, alpha=0.3)

    # Panel D: Inverse selection dynamics
    ax4 = fig.add_subplot(gs[1, 2:4])
    if 'inverse_evolution' in experiment_results:
        gr = experiment_results['inverse_evolution']
        start, step = 0, (gr.generations // 15) + 1
        generations = [i * step for i in range(len(gr.s_means[start::step]))]

        s_means = gr.s_means[start::step]
        q_means = lin_phi_inv(gr.q_means[start::step])
        w_means = gr.w_means[start::step]

        s_norm = (s_means - np.min(s_means)) / (np.max(s_means) - np.min(s_means))
        q_norm = (q_means - np.min(q_means)) / (np.max(q_means) - np.min(q_means))
        w_norm = (w_means - np.min(w_means)) / (np.max(w_means) - np.min(w_means))

        ax4.plot(generations, s_norm, color=colors['blue'], linewidth=3, marker='o',
                markersize=4, label='Structural (s)')
        ax4.plot(generations, q_norm, color=colors['red'], linewidth=3, marker='s',
                markersize=4, label='Quality (q)')
        ax4.plot(generations, w_norm, color=colors['yellow'], linewidth=3, marker='^',
                markersize=4, label='Fitness (w)')

        ax4.set_xlabel('Generation', fontsize=16, fontweight='bold')
        ax4.set_ylabel('Normalized Value', fontsize=16, fontweight='bold')
        ax4.set_title('D. Inverse Selection Dynamics', fontsize=18, fontweight='bold', pad=20)
        ax4.legend(fontsize=14)
        ax4.grid(True, alpha=0.3)

    # Panel E: Neutral evolution
    ax5 = fig.add_subplot(gs[2, 0])
    if 'neutral_evolution' in experiment_results:
        gr = experiment_results['neutral_evolution']
        ax5.text(0.5, 0.5, f'Neutral Evolution\n\nρ(s,q) = {gr.rho_sq_means:.3f}\nQSC = {gr.last_result.QSC:.3f}',
                transform=ax5.transAxes, fontsize=14, ha='center', va='center',
                bbox=dict(boxstyle='round,pad=1', facecolor='lightgray', alpha=0.7))
        ax5.set_title('E. Neutral Selection', fontsize=18, fontweight='bold', pad=20)
        ax5.set_xlim(0, 1)
        ax5.set_ylim(0, 1)
        ax5.axis('off')

    # Panel F: Quality signal limit
    ax6 = fig.add_subplot(gs[2, 1])
    if 'qsl_evolution' in experiment_results:
        gr = experiment_results['qsl_evolution']
        ax6.text(0.5, 0.5, f'Quality Signal Limit\n\nρ(s,q) = {gr.rho_sq_means:.3f}\nQSC = {gr.last_result.QSC:.3f}',
                transform=ax6.transAxes, fontsize=14, ha='center', va='center',
                bbox=dict(boxstyle='round,pad=1', facecolor='lightcoral', alpha=0.7))
        ax6.set_title('F. Signal Limit Scenario', fontsize=18, fontweight='bold', pad=20)
        ax6.set_xlim(0, 1)
        ax6.set_ylim(0, 1)
        ax6.axis('off')

    # Panel G: Mathematical framework summary
    ax7 = fig.add_subplot(gs[2, 2])
    framework_text = r"""Mathematical Framework:

Signal Mapping:
$\phi(s) = 2s^2 + 4 + \epsilon$
$\epsilon \sim \mathcal{N}(0, \sigma^2)$

Selection Dynamics:
$\Delta \bar{z} = \text{Cov}(w, z) + E(w \Delta z)$

Key Parameters:
$\hat{s}$: Signal fidelity (0-1)
QSC: Quality Selection Coefficient
ρ(s,q): Signal correlation"""
    ax7.text(0.05, 0.95, framework_text, transform=ax7.transAxes, fontsize=12,
             verticalalignment='top', fontfamily='monospace',
             bbox=dict(boxstyle='round,pad=0.5', facecolor='lightyellow', alpha=0.8))
    ax7.set_title('G. Mathematical Framework', fontsize=18, fontweight='bold', pad=20)
    ax7.set_xlim(0, 1)
    ax7.set_ylim(0, 1)
    ax7.axis('off')

    # Panel H: Results summary
    ax8 = fig.add_subplot(gs[2, 3])
    results_text = """Key Insights:

• Selection acts on quality signals (q)
  but evolution occurs in structures (s)

• Signal fidelity (\\hat{s}) determines
  selection effectiveness

• Low fidelity leads to signal breakdown
  (Quality Signal Limit)

• Correlated traits can show complex
  evolutionary dynamics (rebound, inverse)

• Fitness landscapes shape signal evolution"""
    ax8.text(0.05, 0.95, results_text, transform=ax8.transAxes, fontsize=12,
             verticalalignment='top',
             bbox=dict(boxstyle='round,pad=0.5', facecolor='lightgreen', alpha=0.8))
    ax8.set_title('H. Key Insights', fontsize=18, fontweight='bold', pad=20)
    ax8.set_xlim(0, 1)
    ax8.set_ylim(0, 1)
    ax8.axis('off')

    # Save the composite figure
    out_path: Path | None = None
    if dest is not None:
        out_path = Path(dest)
        core_io.ensure_directory(out_path.parent)
        plt.savefig(out_path, format='png', dpi=300, bbox_inches='tight',
                   facecolor='white', edgecolor='none')
    plt.close(fig)
    return out_path


