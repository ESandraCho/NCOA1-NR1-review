#!/usr/bin/env python3
"""
gBGC paper — Supplementary figure: NRID GC bias by codon position across clades.

The GC/gBGC bias at the NRID is SYNONYMOUS-confined in the GC-leaning ancestors (reptiles:
elevated at the 3rd/wobble position only, where it does not change the protein), relaxed in
therians, and in monotremes alone OVERFLOWS into the nonsynonymous 1st and 2nd positions — the
positions that determine the amino acid, where it breaks LVQLL -> PGQLP.

This is the same directional bias caught before (reptiles, silent) vs after (monotremes, the
motif lesion) it escapes the synonymous buffer. Complements the main-text branch W->S codon-
position test (Fig 3a).

Reads:  analysis/outgroups_48sp/NCOA1_48sp_codon_aln.fasta (recomputes clade means)
Outputs: figures/FigS_codonpos_gbgc.{png,pdf,svg}
"""
import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams["svg.fonttype"] = "none"
matplotlib.rcParams["font.family"] = "sans-serif"
matplotlib.rcParams["font.sans-serif"] = ["Arial", "Liberation Sans", "Helvetica", "DejaVu Sans"]
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq

PAPER = Path(__file__).resolve().parents[1]
OUT = Path(__file__).resolve().parent

CLADES = {
    "Amphibian": ["Xenopus_tropicalis", "Rana_temporaria", "Bufo_bufo", "Geotrypetes_seraphini"],
    "Reptiles/birds": ["Alligator_mississippiensis", "Crocodylus_porosus", "Anolis_carolinensis",
                       "Gallus_gallus", "Chelonia_mydas", "Podarcis_muralis", "Python_bivittatus",
                       "Sceloporus_undulatus"],
    "Therians": ["Homo_sapiens", "Mus_musculus", "Bos_taurus", "Monodelphis_domestica",
                 "Sarcophilus_harrisii"],
    "Monotremes": ["Ornithorhynchus_anatinus", "Tachyglossus_aculeatus"],
}
COL = {"Amphibian": "#8BC34A", "Reptiles/birds": "#7E57C2", "Therians": "#42A5F5",
       "Monotremes": "#E53935"}


def load():
    return {r.id.replace("_NCOA1", ""): str(r.seq)
            for r in SeqIO.parse(PAPER / "analysis" / "outgroups_48sp" / "NCOA1_48sp_codon_aln.fasta",
                                 "fasta")}


def gc_by_pos(recs, sp, lo, hi):
    s = recs[sp][lo * 3:hi * 3].upper()
    cod = [s[i:i + 3] for i in range(0, len(s), 3) if "-" not in s[i:i + 3]]
    n = len(cod) or 1
    return [100 * sum(c[p] in "GC" for c in cod) / n for p in range(3)]


def main():
    recs = load()
    hum = recs["Homo_sapiens"]
    hp = str(Seq(hum.replace("-", "")).translate())

    def col_of(q):
        mi = hp.find(q); ung = 0
        for c in range(len(hum) // 3):
            if hum[c * 3:c * 3 + 3] != "---":
                if ung == mi:
                    return c
                ung += 1
    nr1 = col_of("LVQLL"); nr3 = col_of("LRYLL")
    lo, hi = nr1 - 12, nr3 + 12

    means = {}
    for cl, sps in CLADES.items():
        vals = [gc_by_pos(recs, s, lo, hi) for s in sps if s in recs]
        means[cl] = [sum(v[p] for v in vals) / len(vals) for p in range(3)]

    fig, ax = plt.subplots(figsize=(6.4, 4.2))
    positions = ["1st\n(nonsyn)", "2nd\n(nonsyn)", "3rd / wobble\n(synonymous)"]
    x = np.arange(3)
    w = 0.2
    for i, (cl, m) in enumerate(means.items()):
        ax.bar(x + (i - 1.5) * w, m, w, label=cl, color=COL[cl], alpha=0.9,
               edgecolor="#333", linewidth=0.4)

    # baseline guide (nonsyn positions sit here in every clade except monotremes)
    ax.axhline(46, color="#999", lw=0.6, ls="--")
    ax.text(2.46, 48, "nonsyn baseline", fontsize=6, color="#777", ha="right", va="bottom")
    # short pointer: reptile bias is synonymous-only (protein intact)
    ax.annotate("reptile bias:\nsynonymous only\n(LVQLL intact)",
                xy=(2 - 1.5 * w + w, means["Reptiles/birds"][2]), xytext=(1.45, 86),
                fontsize=6.5, color="#5E35B1", ha="center",
                arrowprops=dict(arrowstyle="->", color="#5E35B1", lw=0.7))
    # short pointer: monotreme bias invades the nonsynonymous 2nd position
    ax.annotate("monotreme bias\ninvades nonsyn\n→ LVQLL→PGQLP",
                xy=(1 + 1.5 * w, means["Monotremes"][1]), xytext=(0.62, 24),
                fontsize=6.5, color="#E53935", ha="center", fontweight="bold",
                arrowprops=dict(arrowstyle="->", color="#E53935", lw=0.7,
                                connectionstyle="arc3,rad=-0.2"))

    ax.set_xticks(x)
    ax.set_xticklabels(positions, fontsize=8)
    ax.set_ylabel("coding GC content (%)", fontsize=9)
    ax.set_xlabel("Codon position in the NRID", fontsize=9)
    ax.set_title("NRID GC bias is synonymous in ancestors, nonsynonymous only in monotremes",
                 fontsize=9, fontweight="bold")
    ax.set_ylim(0, 100)
    ax.legend(fontsize=7, loc="upper left", framealpha=0.95)
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(labelsize=8)

    for ext in ("png", "pdf", "svg"):
        fig.savefig(OUT / f"FigS_codonpos_gbgc.{ext}", dpi=300, bbox_inches="tight", facecolor="white")
    plt.close()
    print("Saved FigS_codonpos_gbgc")
    for cl, m in means.items():
        print(f"  {cl:16s} pos1={m[0]:.0f} pos2={m[1]:.0f} pos3={m[2]:.0f}")


if __name__ == "__main__":
    main()
