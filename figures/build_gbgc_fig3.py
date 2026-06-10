#!/usr/bin/env python3
"""
gBGC paper — Figure 3: gBGC signature, not adaptation.

(a) W->S vs S->W substitutions on the monotreme stem branch, by codon position. The weak->strong
    bias is confined to nonsynonymous (1st+2nd) sites; synonymous (3rd) sites are not W->S-biased,
    indicating a past, completed conversion rather than ongoing biased fixation, and arguing against
    the elevated-dN/dS being adaptive (the aBSREL/BUSTED p-values are reported in Table S2).
(b) The converted motif is maintained: PGQLP core is 100% identical platypus-echidna; only the
    most GC-extreme flank (NR1) is somewhat less conserved, while NR2/NR3 windows stay ~95%.

Values are computed from data (analysis/gbgc_genomic/compute_reported_values.py), not hard-coded.

Outputs: figures/Fig3_gbgc.{png,pdf,svg}
"""
import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams["svg.fonttype"] = "none"  # emit editable <text>, not glyph paths
matplotlib.rcParams["font.family"] = "sans-serif"
# Arial first (named in SVG; rendered locally via Liberation Sans, the metric-identical clone)
matplotlib.rcParams["font.sans-serif"] = ["Arial", "Liberation Sans", "Helvetica", "DejaVu Sans"]
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import sys
from pathlib import Path

OUT = Path(__file__).resolve().parent
COL_NR1 = "#E53935"

# Pull the plotted numbers from the single canonical computation so the figure and the
# manuscript prose cannot drift apart (previously these were hard-coded in this file).
sys.path.insert(0, str(OUT.parent / "analysis" / "gbgc_genomic"))
from compute_reported_values import ws_counts, platypus_echidna_identities  # noqa: E402


def main():
    fig = plt.figure(figsize=(7.0, 6.0))
    gs = gridspec.GridSpec(2, 1, figure=fig, hspace=0.55, top=0.92, bottom=0.09, left=0.13, right=0.95)

    # (a) W->S by codon position on the stem branch
    ax = fig.add_subplot(gs[0])
    cats = ["Synonymous\n(3rd codon position)", "Nonsynonymous\n(1st + 2nd position)"]
    nws, nsw, sws, ssw = ws_counts()  # computed from codeml ASR, NR1 window
    ws = [sws, nws]
    sw = [ssw, nsw]
    x = np.arange(len(cats)); w = 0.35
    ax.bar(x - w/2, ws, w, label="W→S (A/T→G/C)", color="#5C6BC0")
    ax.bar(x + w/2, sw, w, label="S→W (G/C→A/T)", color="#FF9800")
    for i, (a, b) in enumerate(zip(ws, sw)):
        ax.text(i - w/2, a + 0.15, str(a), ha="center", fontsize=9, fontweight="bold")
        ax.text(i + w/2, b + 0.15, str(b), ha="center", fontsize=9, fontweight="bold")
    ax.set_xticks(x); ax.set_xticklabels(cats, fontsize=9)
    ax.set_ylabel("substitutions on the\nmonotreme stem branch", fontsize=9)
    ax.set_title("Weak→strong bias confined to nonsynonymous sites", fontsize=10, fontweight="bold")
    ax.legend(fontsize=8, loc="upper left")
    ax.text(0, 6.2, "synonymous sites not W→S-biased\n→ a past, completed conversion,\nnot ongoing biased fixation",
            fontsize=7.5, ha="center", color="#666", style="italic")
    ax.set_ylim(0, 12.5)
    ax.spines[["top", "right"]].set_visible(False); ax.tick_params(labelsize=8)

    # (b) retention: motif core vs flank, GC-controlled
    ax = fig.add_subplot(gs[1])
    idd = platypus_echidna_identities()  # computed from the codon alignment
    labels = ["NR1 core\n(PGQLP)", "NR1 flank\n(extreme GC)", "NR2 window\n(elevated GC)",
              "NR3 window\n(elevated GC)", "whole\nprotein"]
    ident = [idd["NR1 core (PGQLP)"], idd["NR1 flank (high GC)"],
             idd["NR2 window (high GC)"], idd["NR3 window (high GC)"], idd["whole protein"]]
    cols = [COL_NR1, "#EF9A9A", "#FFCC80", "#FFCC80", "#90A4AE"]
    ax.bar(range(len(labels)), ident, color=cols, alpha=0.9, edgecolor="#333", linewidth=0.4)
    for i, v in enumerate(ident):
        ax.text(i, v + 1.2, f"{v}%", ha="center", fontsize=9, fontweight="bold")
    ax.set_xticks(range(len(labels))); ax.set_xticklabels(labels, fontsize=8)
    ax.set_ylabel("platypus–echidna identity (%)", fontsize=9)
    ax.set_title("The converted motif is maintained, not drifting", fontsize=10, fontweight="bold")
    ax.set_ylim(0, 112)
    ax.text(2.0, 42, "only the most GC-extreme flank (NR1) is\nless conserved; NR2/NR3 stay ~95% identical",
            fontsize=7.5, ha="center", color="#666", style="italic")
    ax.spines[["top", "right"]].set_visible(False); ax.tick_params(labelsize=8)

    fig.text(0.02, 0.93, "(a)", fontsize=14, fontweight="bold")
    fig.text(0.02, 0.46, "(b)", fontsize=14, fontweight="bold")
    for ext in ("png", "pdf", "svg"):
        fig.savefig(OUT / f"Fig3_gbgc.{ext}", dpi=300, bbox_inches="tight", facecolor="white")
    plt.close()
    print("Saved Fig3_gbgc (2-panel: W->S bias + retention)")


if __name__ == "__main__":
    main()
