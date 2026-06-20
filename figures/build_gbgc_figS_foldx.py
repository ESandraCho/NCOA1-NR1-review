#!/usr/bin/env python3
"""
gBGC paper — Supplementary Figure S3: the overrun NR1 residue is structurally consequential.

FoldX (AnalyseComplex) interaction energy of the NCOA1 NR-box peptide bound to its receptor
(co-crystal PDB 1XJ7), for three peptide sequences modelled in place (BuildModel, 5 models each):
  - LVQLL  : ancestral / non-monotreme NR1 (the conserved state)
  - LHRLL  : intact NR2 sequence, included as a control for a different but functional NR box
  - PGQLP  : the monotreme NR1 state produced by gBGC

Result: ancestral LVQLL and control LHRLL bind strongly (~ -15 kcal/mol), while the monotreme
PGQLP binds markedly more weakly (~ -8 kcal/mol; ΔΔG ≈ +7 kcal/mol). The two introduced prolines
(helix-breaking) disrupt the amphipathic docking helix. This is a STRUCTURAL PREDICTION, not a
functional measurement — FoldX overestimates proline effects — and is interpreted only as evidence
that the ancestral residue position contributed substantially to docking (i.e. gBGC overran a
consequential, not a tolerant, position).

Reads:   analysis/foldx/Interaction_1XJ7_Repair_{LVQLL,LHRLL,PGQLP}_*_AC.fxout
Outputs: figures/FigS_foldx_gbgc.{png,pdf,svg}
"""
import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams["svg.fonttype"] = "none"
matplotlib.rcParams["font.family"] = "sans-serif"
matplotlib.rcParams["font.sans-serif"] = ["Arial", "Liberation Sans", "Helvetica", "DejaVu Sans"]
import matplotlib.pyplot as plt
import numpy as np
import glob
from pathlib import Path

PAPER = Path(__file__).resolve().parents[1]
FOLDX = PAPER / "analysis" / "foldx"
OUT = Path(__file__).resolve().parent

COL_NR1 = "#E53935"   # red — the disrupted monotreme state
COL_ANC = "#90A4AE"   # grey — ancestral / intact

# (sequence key, display label, bar colour)
VARIANTS = [
    ("LVQLL", "LVQLL\n(ancestral NR1)", COL_ANC),
    ("LHRLL", "LHRLL\n(NR2 control)", COL_ANC),
    ("PGQLP", "PGQLP\n(monotreme NR1)", COL_NR1),
]


def interaction_energy(path):
    """Interaction Energy (column 6) from a FoldX AnalyseComplex .fxout file."""
    for ln in Path(path).read_text().splitlines():
        cols = ln.split("\t")
        if cols and cols[0].endswith(".pdb") and len(cols) > 5:
            try:
                return float(cols[5])
            except ValueError:
                pass
    return None


def variant_values(key):
    files = sorted(glob.glob(str(FOLDX / f"Interaction_1XJ7_Repair_{key}_*_AC.fxout")))
    vals = [interaction_energy(f) for f in files]
    return [v for v in vals if v is not None]


def main():
    data = {k: variant_values(k) for k, _, _ in VARIANTS}
    means = {k: float(np.mean(v)) for k, v in data.items()}
    sds = {k: float(np.std(v)) for k, v in data.items()}
    ddg = means["PGQLP"] - means["LVQLL"]

    fig, ax = plt.subplots(figsize=(4.8, 4.0))
    x = np.arange(len(VARIANTS))
    labels = [lab for _, lab, _ in VARIANTS]
    cols = [c for _, _, c in VARIANTS]
    mvals = [means[k] for k, _, _ in VARIANTS]
    evals = [sds[k] for k, _, _ in VARIANTS]

    ax.bar(x, mvals, color=cols, alpha=0.9, edgecolor="#333", linewidth=0.5,
           yerr=evals, capsize=4, error_kw=dict(lw=0.8, ecolor="#333"))
    # overlay the 5 per-model points per variant
    for i, (k, _, _) in enumerate(VARIANTS):
        jitter = (np.arange(len(data[k])) - (len(data[k]) - 1) / 2) * 0.05
        ax.scatter(np.full(len(data[k]), x[i]) + jitter, data[k],
                   s=12, color="#222", zorder=3, alpha=0.7)
    # value labels, placed just below each bar top (inside the bar, dark text)
    for i, (k, _, _) in enumerate(VARIANTS):
        ax.text(x[i] + 0.30, means[k], f"{means[k]:.1f}", ha="left", va="center",
                fontsize=8, fontweight="bold", color="#222")

    # ΔΔG annotation: a horizontal level marker from the ancestral bar to the monotreme bar top
    ax.annotate("", xy=(2.0, means["PGQLP"]), xytext=(0.0, means["PGQLP"]),
                arrowprops=dict(arrowstyle="-", color="#888", lw=0.7, ls=":"))
    ax.annotate("", xy=(2.0, means["LVQLL"]), xytext=(2.0, means["PGQLP"]),
                arrowprops=dict(arrowstyle="<->", color="#555", lw=0.9))
    ax.text(2.12, (means["LVQLL"] + means["PGQLP"]) / 2,
            f"ΔΔG ≈ +{ddg:.1f}\nkcal/mol", ha="left", va="center", fontsize=7.5,
            color="#555", style="italic")

    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=8)
    ax.set_ylabel("FoldX interaction energy (kcal/mol)\nmore negative = stronger binding", fontsize=9)
    ax.set_title("Monotreme NR1 (PGQLP) is predicted to weaken receptor docking",
                 fontsize=8.2, fontweight="bold")
    ax.set_ylim(min(mvals) - 2.5, 0)
    ax.set_xlim(-0.6, 2.95)
    ax.axhline(0, color="#999", lw=0.6)
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(labelsize=8)
    fig.tight_layout()

    for ext in ("png", "pdf", "svg"):
        fig.savefig(OUT / f"FigS_foldx_gbgc.{ext}", dpi=300, bbox_inches="tight", facecolor="white")
    plt.close()

    print("Saved FigS_foldx_gbgc")
    for k, _, _ in VARIANTS:
        print(f"  {k}: mean {means[k]:.2f} kcal/mol  (n={len(data[k])}, sd={sds[k]:.2f})")
    print(f"  ΔΔG (PGQLP - LVQLL) = {ddg:+.2f} kcal/mol")


if __name__ == "__main__":
    main()
