#!/usr/bin/env python3
"""
gBGC paper — VERIFICATION (not a figure): the PGQLP-for-LVQLL replacement is alignment-robust.

This result is reported as text in the manuscript (Methods/Results), not as a supplementary figure
— the informational content is a single statement ("the call is identical across three alignment
methods"), so it does not warrant a panel. This script is the reproducible check behind that claim:
it re-aligns the NCOA1 proteins under three independent MAFFT regimes — L-INS-i
(--localpair --maxiterate 1000), the FFT-NS-2 default (--retree 2), and a high gap-open penalty
(--op 3.0) — and prints the NR1 column for human and the two monotremes. The motif is LVQLL (human)
and PGQLP (both monotremes) in all three regimes, confirming the disruption is not an alignment
artifact. (It can still render a box panel if --plot is passed, but no rendered figure is committed.)

Reads:  analysis/outgroups_48sp/NCOA1_48sp_prot.fasta (re-aligns under the 3 regimes)
Requires: mafft on PATH.
"""
import sys
import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams["svg.fonttype"] = "none"
matplotlib.rcParams["font.family"] = "sans-serif"
matplotlib.rcParams["font.sans-serif"] = ["Arial", "Liberation Sans", "Helvetica", "DejaVu Sans"]
import matplotlib.pyplot as plt
import subprocess
import tempfile
import os
from pathlib import Path
from Bio import SeqIO

PAPER = Path(__file__).resolve().parents[1]
OUT = Path(__file__).resolve().parent
COL_NR1 = "#E53935"
COL_OK = "#2E7D32"
BG_OK = "#C8E6C9"
BG_BAD = "#FFCDD2"

REGIMES = {
    "MAFFT L-INS-i": ["mafft", "--localpair", "--maxiterate", "1000", "--quiet"],
    "MAFFT FFT-NS-2 (default)": ["mafft", "--retree", "2", "--quiet"],
    "MAFFT high gap penalty (op=3)": ["mafft", "--op", "3.0", "--quiet"],
}
SPECIES = [("Human (therian)", "Homo_sapiens"),
           ("Platypus", "Ornithorhynchus_anatinus"),
           ("Echidna", "Tachyglossus_aculeatus")]


def align(cmd, infile):
    out = subprocess.run(cmd + [infile], capture_output=True, text=True, check=True).stdout
    a, n = {}, None
    for ln in out.splitlines():
        if ln.startswith(">"):
            n = ln[1:].split()[0]; a[n] = ""
        elif n:
            a[n] += ln.strip()
    return a


def nr1_col(aln):
    h = aln["Homo_sapiens"]; ung = h.replace("-", ""); mi = ung.find("LVQLL"); seen = 0
    for c, ch in enumerate(h):
        if ch != "-":
            if seen == mi:
                return c
            seen += 1


def main():
    prot = {r.id.replace("_NCOA1", ""): str(r.seq).replace("-", "")
            for r in SeqIO.parse(PAPER / "analysis" / "outgroups_48sp" / "NCOA1_48sp_prot.fasta", "fasta")}
    with tempfile.NamedTemporaryFile("w", suffix=".fa", delete=False) as f:
        for k, v in prot.items():
            f.write(f">{k}\n{v}\n")
        inp = f.name

    # collect NR1 5-mer per (regime, species)
    table = {}
    for rname, cmd in REGIMES.items():
        aln = align(cmd, inp)
        col = nr1_col(aln)
        table[rname] = {sp: aln[sp][col:col + 5] for _, sp in SPECIES}
    os.unlink(inp)

    # Always print the verification result (this is the point of the script).
    print("Alignment-robustness check — NR1 column per MAFFT regime:")
    for rname in REGIMES:
        print(f"  {rname:32s} " + "  ".join(f"{lbl.split()[0]}={table[rname][s]}"
                                            for lbl, s in SPECIES))
    if "--plot" not in sys.argv:
        print("(no figure committed; pass --plot to render an optional box panel)")
        return

    fig, ax = plt.subplots(figsize=(8.6, 3.4))
    ax.set_xlim(0, 16); ax.set_ylim(0, len(REGIMES) + 1.3); ax.axis("off")
    ax.text(8, len(REGIMES) + 0.9, "NR1 column under three alignment regimes",
            fontsize=10, ha="center", fontweight="bold")

    cw = 0.62
    # left margin 4.6 for the regime labels; each species block is 5*cw wide (~3.1) + a gap
    x0 = {sp_label: 4.8 + i * 3.6 for i, (sp_label, _) in enumerate(SPECIES)}
    # species headers
    for sp_label, _ in SPECIES:
        ax.text(x0[sp_label] + 2 * cw, len(REGIMES) + 0.25, sp_label, fontsize=8,
                ha="center", color="#555", fontweight="bold")

    for r, (rname, _) in enumerate(REGIMES.items()):
        y = len(REGIMES) - r - 0.3
        ax.text(0.1, y, rname, fontsize=8, va="center")
        for sp_label, sp in SPECIES:
            mot = table[rname][sp]
            dis = mot != "LVQLL"
            for j, aa in enumerate(mot):
                x = x0[sp_label] + j * cw
                ax.add_patch(plt.Rectangle((x - cw / 2 + 0.03, y - 0.32), cw - 0.06, 0.64,
                                           facecolor=BG_BAD if dis else BG_OK, edgecolor="none",
                                           alpha=0.75))
                ax.text(x, y, aa, fontsize=8.5, family="monospace", ha="center", va="center",
                        color=COL_NR1 if dis else COL_OK, fontweight="bold")

    ax.text(8, 0.15, "Human LVQLL and monotreme PGQLP are identical across all three regimes — "
                     "the disruption is not an alignment artifact.",
            fontsize=6.5, ha="center", style="italic", color="#444")

    for ext in ("png", "pdf", "svg"):
        fig.savefig(OUT / f"FigS_alnrobust_gbgc.{ext}", dpi=300, bbox_inches="tight", facecolor="white")
    plt.close()
    print("Rendered optional panel FigS_alnrobust_gbgc.{png,pdf,svg} (not committed).")


if __name__ == "__main__":
    main()
