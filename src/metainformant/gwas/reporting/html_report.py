"""Self-contained HTML report generator for GWAS post-analysis.

Produces a single-file HTML summary with inline CSS — no external
dependencies required for rendering.  Embeds statistics cards,
gene annotation tables, visualization gallery, and supplementary
tests (sign test, credible sets, conditional analysis).
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List, Optional

from metainformant.core.utils.logging import get_logger

logger = get_logger(__name__)


def _fmt_p(p: float) -> str:
    """Format a p-value for display."""
    try:
        return f"{float(p):.3e}"
    except Exception:
        return str(p)


def generate_html_report(
    output_dir: str | Path,
    summary_stats: Dict[str, Any],
    *,
    sign_test_results: Optional[Dict[str, Any]] = None,
    credible_set_results: Optional[Dict[str, Any]] = None,
    conditional_results: Optional[Dict[str, Any]] = None,
    gene_annotations: Optional[List[Dict[str, Any]]] = None,
    trait: str = "unknown",
    model: str = "unknown",
    n_samples: Optional[int] = None,
    filename: str = "report.html",
) -> Path:
    """Generate a self-contained HTML summary report.

    Args:
        output_dir: Directory to write the report into.
        summary_stats: Dict from ``compute_comprehensive_summary``.
        sign_test_results: Optional sign-test dict.
        credible_set_results: Optional credible-set dict.
        conditional_results: Optional conditional analysis dict.
        gene_annotations: Optional list of gene annotations.
        trait: Trait name for the report header.
        model: GWAS model name.
        n_samples: Number of samples.
        filename: Output filename (default ``report.html``).

    Returns:
        Path to the generated HTML file.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    ss = summary_stats
    sig = ss.get("significance_counts", {})
    maf = ss.get("maf_distribution", {})
    eff = ss.get("effect_size_distribution", {})
    cal = ss.get("p_value_calibration", {})
    top_hits = ss.get("top_10_hits", [])

    # ── Sign test row ──
    sign_row = ""
    if sign_test_results:
        st = sign_test_results
        sign_row = f"""
        <tr><td>Effect direction (sign test)</td>
        <td>n+={st.get("n_positive_beta", "?")}, n−={st.get("n_negative_beta", "?")},
        z={st.get("z_statistic", "?")}, p={_fmt_p(st.get("p_value_two_sided", 1))}</td>
        <td>{"✅" if st.get("p_value_two_sided", 1) < 0.05 else "—"}</td></tr>
        <tr><td colspan="3" style="color:#555;font-size:0.9em;padding-left:20px">
        {st.get("interpretation", "")}</td></tr>
        """

    # ── Credible set row ──
    cs_row = ""
    if credible_set_results:
        cs = credible_set_results
        cs_row = f"""<tr><td>Credible set (95% coverage)</td>
        <td>{cs.get("n_snps", "?")} variants</td><td>✅</td></tr>"""

    # ── Conditional analysis row ──
    cond_row = ""
    if conditional_results:
        cr = conditional_results
        cond_row = f"""<tr><td>Independent signals (conditional)</td>
        <td>{cr.get("n_independent_signals", "?")}</td><td>✅</td></tr>"""

    # ── Gene annotation rows ──
    gene_rows = ""
    for ann in (gene_annotations or [])[:5]:
        genes = ann.get("nearby_genes", [])
        gene_str = (
            ", ".join([g.get("gene_name", "?") for g in genes[:3]]) if genes else "—"
        )
        gene_rows += f"""<tr>
            <td>{ann.get("snp", "?")}</td>
            <td>{ann.get("chrom", "?")}:{ann.get("pos", "?")}</td>
            <td>{gene_str}</td>
        </tr>"""

    # ── Top hits rows ──
    top_hit_rows = ""
    for h in top_hits[:5]:
        top_hit_rows += f"""<tr>
            <td>{h.get("snp", "?")}</td>
            <td>{h.get("chrom", "?")}:{h.get("pos", "?"):,}</td>
            <td>{h.get("beta", "?"):.3f}</td>
            <td>{h.get("se", "?"):.3f}</td>
            <td>{_fmt_p(h.get("p_value", 1))}</td>
            <td>{h.get("MAF", "?"):.3f}</td>
        </tr>"""

    lambda_gc = cal.get("lambda_gc", float("nan"))
    lambda_warn = (
        "⚠️ Inflated λ_GC — check for stratification or strong signal"
        if lambda_gc > 2
        else ""
    )

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>GWAS Report — {trait} ({model})</title>
<style>
  body {{ font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Helvetica, sans-serif;
         max-width: 1100px; margin: 40px auto; color: #1a1a1a; line-height: 1.6; }}
  h1 {{ font-size: 2em; border-bottom: 3px solid #0057a8; padding-bottom: 8px; color: #003f7a; }}
  h2 {{ font-size: 1.3em; color: #0057a8; margin-top: 2em; }}
  table {{ border-collapse: collapse; width: 100%; margin: 1em 0; font-size: 0.92em; }}
  th {{ background: #0057a8; color: white; padding: 8px 12px; text-align: left; }}
  td {{ padding: 6px 12px; border-bottom: 1px solid #ddd; }}
  tr:hover {{ background: #f0f4ff; }}
  .stat-grid {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(160px, 1fr));
               gap: 12px; margin: 1em 0; }}
  .stat-card {{ background: linear-gradient(135deg, #f0f7ff, #e3eeff);
                border-radius: 8px; padding: 14px 16px; border-left: 4px solid #0057a8; }}
  .stat-card .value {{ font-size: 1.6em; font-weight: 700; color: #003f7a; }}
  .stat-card .label {{ font-size: 0.82em; color: #555; margin-top: 2px; }}
  .warn {{ background: #fff3cd; border-left: 4px solid #ffc107; padding: 8px 14px;
           border-radius: 4px; font-size: 0.9em; margin: 8px 0; }}
  .plot-grid {{ display: grid; grid-template-columns: 1fr 1fr; gap: 20px; }}
  .plot-grid img {{ width: 100%; border-radius: 6px; border: 1px solid #ddd; }}
  footer {{ margin-top: 40px; font-size: 0.8em; color: #888; border-top: 1px solid #eee;
            padding-top: 12px; }}
</style>
</head>
<body>
<h1>GWAS Summary Report</h1>
<p><strong>Trait:</strong> {trait} &nbsp;|&nbsp;
   <strong>Model:</strong> {model} &nbsp;|&nbsp;
   <strong>Variants tested:</strong> {ss.get("n_variants", "?"):,} &nbsp;|&nbsp;
   <strong>Chromosomes:</strong> {ss.get("n_chromosomes", "?")}</p>

<h2>Key Statistics</h2>
<div class="stat-grid">
  <div class="stat-card">
    <div class="value">{sig.get("genome_wide_5e-8", 0)}</div>
    <div class="label">GW-significant hits (p&lt;5e-8)</div>
  </div>
  <div class="stat-card">
    <div class="value">{sig.get("suggestive_1e-5", 0)}</div>
    <div class="label">Suggestive hits (p&lt;1e-5)</div>
  </div>
  <div class="stat-card">
    <div class="value">{lambda_gc:.3f}</div>
    <div class="label">λ_GC (genomic inflation)</div>
  </div>
  <div class="stat-card">
    <div class="value">{_fmt_p(cal.get("min_p", 1))}</div>
    <div class="label">Minimum p-value</div>
  </div>
  <div class="stat-card">
    <div class="value">{maf.get("mean", 0):.3f}</div>
    <div class="label">Mean MAF</div>
  </div>
  <div class="stat-card">
    <div class="value">{eff.get("mean_abs_beta", 0):.3f}</div>
    <div class="label">Mean |β|</div>
  </div>
</div>
{'<div class="warn">' + lambda_warn + "</div>" if lambda_warn else ""}

<h2>Top 5 Association Hits</h2>
<table>
  <tr><th>SNP</th><th>Location</th><th>β</th><th>SE</th><th>P-value</th><th>MAF</th></tr>
  {top_hit_rows}
</table>

<h2>Supplementary Statistics</h2>
<table>
  <tr><th>Test</th><th>Result</th><th>Pass</th></tr>
  {sign_row}
  {cs_row}
  {cond_row}
  <tr><td>Calibration (q50 obs/exp ratio)</td>
      <td>{cal.get("calibration_quantiles", dict()).get("q50", dict()).get("ratio_obs_exp", "—")}</td>
      <td>{"✅" if abs(cal.get("calibration_quantiles", dict()).get("q50", dict()).get("ratio_obs_exp", 1) - 1) < 0.5 else "⚠️"}</td></tr>
</table>

{"<h2>Gene Proximity Annotation (Top 5 Hits)</h2><table><tr><th>SNP</th><th>Location</th><th>Nearby Genes</th></tr>" + gene_rows + "</table>" if gene_rows else ""}

<h2>Visualization Gallery</h2>
<div class="plot-grid">
  <div><p><strong>Manhattan Plot</strong></p><img src="../../manhattan_plot.png" alt="Manhattan plot"></div>
  <div><p><strong>Q-Q Plot</strong></p><img src="../../qq_plot.png" alt="QQ plot"></div>
  <div><p><strong>Volcano Plot</strong></p><img src="volcano_plot.png" alt="Volcano plot"></div>
  <div><p><strong>Stratified Q-Q</strong></p><img src="qq_stratified.png" alt="Stratified QQ"></div>
  <div><p><strong>LD Heatmap (Top Locus)</strong></p><img src="ld_heatmap.png" alt="LD heatmap"></div>
  <div><p><strong>Missingness</strong></p><img src="missingness.png" alt="Missingness"></div>
  <div><p><strong>PCA (Population Structure)</strong></p><img src="../../pca_plot.png" alt="PCA"></div>
  <div><p><strong>Kinship Matrix</strong></p><img src="../../kinship_plot.png" alt="Kinship"></div>
  <div><p><strong>Heritability (per chromosome)</strong></p><img src="../../heritability_bar_chart.png" alt="Heritability"></div>
  <div><p><strong>Effect Size Distribution</strong></p><img src="../../effect_size_plot.png" alt="Effect sizes"></div>
  <div><p><strong>MAF Spectrum</strong></p><img src="../../maf_spectrum_plot.png" alt="MAF spectrum"></div>
  <div><p><strong>Power Curves</strong></p><img src="../power_analysis/power_curves.png" alt="Power curves"></div>
</div>

<footer>Generated by MetaInformAnt GWAS Pipeline &bull; Trait: {trait} &bull; Model: {model}</footer>
</body>
</html>
"""
    report_path = output_dir / filename
    with open(report_path, "w", encoding="utf-8") as f:
        f.write(html)
    logger.info("HTML report saved to %s", report_path)
    return report_path
