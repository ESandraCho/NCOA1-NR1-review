#!/usr/bin/env python3
"""
gBGC paper — Supplementary Figure S1: the NRID CpG/GC enrichment is graded, monotreme-extreme.

Tests whether the NCOA1 NRID CpG enrichment is monotreme-unique or part of a broader gradient. For
each species' genomic NCOA1 locus, the observed/expected CpG ratio (Gardiner-Garden) of the ±750 bp
NR1 window is compared with the local genomic background.

Result: the enrichment is GRADED. GC-leaning sauropsids are already elevated above their own
background (anole 0.52, chicken 0.53), the monotremes higher still (echidna 0.57; platypus 0.74 —
the only locus clearing the formal o/e>0.6 island threshold), while the amphibian, marsupials and
placental sit at the depleted baseline (0.15–0.34). The NRID enrichment is thus an ancestral
GC-leaning feature intensified to the island threshold in the monotreme lineage, NOT a therian
loss — consistent with the "highest in monotremes" GC framing and a monotreme gBGC/recombination
intensification on a pre-existing substrate.

Reads:  analysis/gbgc_genomic/{species}_NCOA1_genomic.fasta
Outputs: figures/FigS_cpgisland_gbgc.{png,pdf,svg}
"""
import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams["svg.fonttype"] = "none"
matplotlib.rcParams["font.family"] = "sans-serif"
matplotlib.rcParams["font.sans-serif"] = ["Arial", "Liberation Sans", "Helvetica", "DejaVu Sans"]
import matplotlib.pyplot as plt
import numpy as np
import sys
from pathlib import Path
from Bio import SeqIO

PAPER = Path(__file__).resolve().parents[1]
OUT = Path(__file__).resolve().parent
COL_NR1 = "#E53935"

# Locate the NR1 exon by motif anchoring (handles both strands; not a hard-coded offset).
sys.path.insert(0, str(PAPER / "analysis" / "gbgc_genomic"))
from compute_reported_values import orient_and_find  # noqa: E402

# (display label, file key, clade colour).
PANEL = [
    ("frog (amphibian)", "frog", "#8BC34A"),
    ("anole (reptile)", "anole", "#7E57C2"),
    ("turtle (reptile)", "turtle", "#7E57C2"),
    ("crocodile (reptile)", "crocodile", "#7E57C2"),
    ("chicken (bird)", "chicken", "#7E57C2"),
    ("platypus (monotreme)", "platypus", COL_NR1),
    ("echidna (monotreme)", "echidna", COL_NR1),
    ("opossum (marsupial)", "opossum", "#26A69A"),
    ("devil (marsupial)", "devil", "#26A69A"),
    ("human (placental)", "human", "#42A5F5"),
]
def oe_cpg(s):
    s = s.upper(); n = len(s)
    if n < 2:
        return 0.0
    cpg = s.count("CG"); C = s.count("C"); G = s.count("G")
    exp = (C * G) / n
    return cpg / exp if exp > 0 else 0.0


def load_oriented(key):
    raw = str(SeqIO.read(PAPER / "analysis" / "gbgc_genomic" / f"{key}_NCOA1_genomic.fasta",
                         "fasta").seq).upper()
    seq, off = orient_and_find(raw)
    if off is None:
        raise RuntimeError(f"NR1 exon not located in {key} genomic locus")
    return seq, off


def main():
    labels, nrid, bg, cols = [], [], [], []
    for label, key, col in PANEL:
        g, pos = load_oriented(key)
        nrid_seq = g[pos - 750:pos + 750]
        bg_seq = g[:pos - 1500] + g[pos + 1500:]
        labels.append(label); nrid.append(oe_cpg(nrid_seq)); bg.append(oe_cpg(bg_seq)); cols.append(col)

    x = np.arange(len(labels))
    fig, ax = plt.subplots(figsize=(8.4, 4.4))
    w = 0.38
    ax.bar(x - w / 2, nrid, w, color=cols, alpha=0.9, edgecolor="#333", linewidth=0.4,
           label="NRID window (±750 bp)")
    ax.bar(x + w / 2, bg, w, color="#CCCCCC", alpha=0.85, edgecolor="#333", linewidth=0.4,
           label="local genomic background")
    # CpG-island threshold
    ax.axhline(0.6, color="#E53935", lw=0.8, ls="--")
    ax.text(len(labels) - 0.5, 0.62, "CpG-island threshold (o/e > 0.6)", fontsize=6.5,
            color="#E53935", ha="right", va="bottom")
    # mark which clear the island criterion
    for i, (label, key, col) in enumerate(PANEL):
        if nrid[i] > 0.6:
            ax.text(i - w / 2, nrid[i] + 0.02, "island", fontsize=6, ha="center",
                    color=COL_NR1, fontweight="bold")

    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=7, rotation=35, ha="right")
    ax.set_ylabel("observed / expected CpG", fontsize=9)
    ax.set_title("NRID CpG/GC enrichment is graded across amniotes, monotreme-extreme (only platypus clears 0.6)",
                 fontsize=7.8, fontweight="bold")
    ax.set_ylim(0, 0.9)
    ax.legend(fontsize=7, loc="upper left", framealpha=0.95)
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(labelsize=7)
    fig.tight_layout()

    for ext in ("png", "pdf", "svg"):
        fig.savefig(OUT / f"FigS_cpgisland_gbgc.{ext}", dpi=300, bbox_inches="tight", facecolor="white")
    plt.close()
    print("Saved FigS_cpgisland_gbgc")
    for label, n, b in zip(labels, nrid, bg):
        print(f"  {label:22s} NRID o/e={n:.2f}  bkg={b:.2f}  {'ISLAND' if n>0.6 else ''}")


if __name__ == "__main__":
    main()
