#!/usr/bin/env python3
"""
gBGC paper — Figure 1: a 400-My-conserved coactivator motif uniquely lost in monotremes.

(a) NR1 motif (LXXLL) across 48 sarcopterygians by clade: LVQLL invariant from coelacanth/
    amphibians through reptiles, birds and therians; PGQLP only in the two monotremes. Marker
    at the monotreme node shows the ancestral reconstruction (PGQLP, AA posterior ≥0.999).
(b) The three NR boxes: NR2 (LHRLL) and NR3 (LRYLL) are intact in monotremes; only NR1 changed.

Reads: analysis/outgroups_48sp/NCOA1_48sp_codon_aln.fasta
Outputs: figures/Fig1_gbgc.{png,pdf,svg}
"""
import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams["svg.fonttype"] = "none"  # emit editable <text>, not glyph paths
matplotlib.rcParams["font.family"] = "sans-serif"
# Arial first (named in SVG; rendered locally via Liberation Sans, the metric-identical clone)
matplotlib.rcParams["font.sans-serif"] = ["Arial", "Liberation Sans", "Helvetica", "DejaVu Sans"]
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import sys
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq

PAPER = Path(__file__).resolve().parents[1]
OUT = Path(__file__).resolve().parent

# Ancestral NR1 motif at the monotreme stem, read from the codeml rst (not typed).
sys.path.insert(0, str(PAPER / "analysis" / "gbgc_genomic"))
from compute_reported_values import stem_nr1_motif  # noqa: E402
COL_NR1 = "#E53935"
COL_OK = "#2E7D32"
BG_OK = "#C8E6C9"
BG_BAD = "#FFCDD2"

# representative species per clade (subset for legibility), in display order
CLADES = [
    ("Outgroups", ["Latimeria_chalumnae", "Geotrypetes_seraphini", "Bufo_bufo", "Xenopus_tropicalis"]),
    ("Reptiles/birds", ["Chelonia_mydas", "Gallus_gallus", "Crocodylus_porosus", "Anolis_carolinensis",
                        "Python_bivittatus"]),
    ("Monotremes", ["Ornithorhynchus_anatinus", "Tachyglossus_aculeatus"]),
    ("Marsupials", ["Monodelphis_domestica", "Sarcophilus_harrisii"]),
    ("Eutherians", ["Homo_sapiens", "Mus_musculus", "Bos_taurus"]),
]
COMMON = {"Latimeria_chalumnae": "coelacanth", "Geotrypetes_seraphini": "caecilian",
          "Bufo_bufo": "toad", "Xenopus_tropicalis": "frog", "Chelonia_mydas": "turtle",
          "Gallus_gallus": "chicken", "Crocodylus_porosus": "crocodile", "Anolis_carolinensis": "anole",
          "Python_bivittatus": "python", "Ornithorhynchus_anatinus": "platypus",
          "Tachyglossus_aculeatus": "echidna", "Monodelphis_domestica": "opossum",
          "Sarcophilus_harrisii": "devil", "Homo_sapiens": "human", "Mus_musculus": "mouse",
          "Bos_taurus": "cow"}

FULL = {'TTT':'F','TTC':'F','TTA':'L','TTG':'L','CTT':'L','CTC':'L','CTA':'L','CTG':'L','ATT':'I','ATC':'I','ATA':'I','ATG':'M','GTT':'V','GTC':'V','GTA':'V','GTG':'V','TCT':'S','TCC':'S','TCA':'S','TCG':'S','CCT':'P','CCC':'P','CCA':'P','CCG':'P','ACT':'T','ACC':'T','ACA':'T','ACG':'T','GCT':'A','GCC':'A','GCA':'A','GCG':'A','TAT':'Y','TAC':'Y','CAT':'H','CAC':'H','CAA':'Q','CAG':'Q','AAT':'N','AAC':'N','AAA':'K','AAG':'K','GAT':'D','GAC':'D','GAA':'E','GAG':'E','TGT':'C','TGC':'C','TGG':'W','CGT':'R','CGC':'R','CGA':'R','CGG':'R','AGT':'S','AGC':'S','AGA':'R','AGG':'R','GGT':'G','GGC':'G','GGA':'G','GGG':'G','---':'-'}


def load():
    return {r.id.replace("_NCOA1", ""): str(r.seq)
            for r in SeqIO.parse(PAPER / "analysis" / "outgroups_48sp" / "NCOA1_48sp_codon_aln.fasta", "fasta")}


def motif_at(recs, motif_aa, sp_seq_col):
    pass


def find_box_col(recs, query):
    hum = recs["Homo_sapiens"]
    hp = str(Seq(hum.replace("-", "")).translate())
    mi = hp.find(query)
    ung, col = 0, None
    for c in range(len(hum) // 3):
        if hum[c*3:c*3+3] != "---":
            if ung == mi:
                col = c; break
            ung += 1
    return col


def motif(recs, sp, col):
    return "".join(FULL.get(recs[sp][i*3:i*3+3].upper(), 'x') for i in range(col, col+5))


def panel_a(ax, recs):
    nr1 = find_box_col(recs, "LVQLL")
    # flatten species in display order with clade headers
    rows = []
    for clade, sps in CLADES:
        rows.append(("HEADER", clade, None))
        for sp in sps:
            if sp in recs:
                rows.append(("SP", sp, motif(recs, sp, nr1)))
    n = len(rows)
    ax.set_xlim(0, 12); ax.set_ylim(-1, n + 1); ax.axis("off")
    ax.text(6, n + 0.3, "NR box 1", fontsize=9, ha="center", fontweight="bold")
    cw = 0.7; mx = 6.5
    for i, (kind, name, mot) in enumerate(rows):
        y = n - 1 - i
        if kind == "HEADER":
            ax.text(0.1, y, name, fontsize=8, fontweight="bold", color="#444", va="center")
            continue
        disp = name.replace("_", " ")
        gs = disp.split()
        lab = f"{gs[0][0]}. {' '.join(gs[1:])}"
        col = COL_NR1 if mot != "LVQLL" else "#666"
        ax.text(0.6, y, lab, fontsize=6, fontstyle="italic", va="center", color=col)
        ax.text(3.6, y, COMMON.get(name, ""), fontsize=5, va="center", color="#aaa")
        dis = mot != "LVQLL"
        for j, aa in enumerate(mot):
            x = mx + j * cw
            ax.add_patch(plt.Rectangle((x - cw/2 + 0.03, y - 0.38), cw - 0.06, 0.76,
                                       facecolor=BG_BAD if dis else BG_OK, edgecolor="none", alpha=0.7))
            ax.text(x, y, aa, fontsize=7.5, family="monospace", ha="center", va="center",
                    color=COL_NR1 if dis else COL_OK, fontweight="bold")
    # ASR annotation at monotreme block
    mono_ys = [n - 1 - i for i, (k, nm, m) in enumerate(rows) if k == "SP" and nm in
               ("Ornithorhynchus_anatinus", "Tachyglossus_aculeatus")]
    if mono_ys:
        ymid = sum(mono_ys) / len(mono_ys)
        ax.annotate(f"ancestral NR1\n{stem_nr1_motif()} (PP ≥ 0.999)", xy=(mx - 0.5, ymid),
                    xytext=(1.5, ymid + 2.2), fontsize=5.5, ha="center", color=COL_NR1, fontweight="bold",
                    arrowprops=dict(arrowstyle="-", color=COL_NR1, lw=0.6))
    ax.set_title("NR1 (LVQLL) invariant across 46 species; PGQLP only in monotremes",
                 fontsize=8.5, fontweight="bold", pad=14)


def panel_b(ax, recs):
    # three NR boxes for human vs the two monotremes
    boxes = [("NR1", find_box_col(recs, "LVQLL")),
             ("NR2", find_box_col(recs, "LHRLL")),
             ("NR3", find_box_col(recs, "LRYLL"))]
    sps = [("Therian (human)", "Homo_sapiens"), ("Platypus", "Ornithorhynchus_anatinus"),
           ("Echidna", "Tachyglossus_aculeatus")]
    ax.set_xlim(0, 12); ax.set_ylim(0, 5); ax.axis("off")
    ax.text(6, 4.6, "Only NR1 changed; NR2/NR3 intact in monotremes", fontsize=8.5,
            ha="center", fontweight="bold")
    cw = 0.62
    for r, (lab, sp) in enumerate(sps):
        y = 3.4 - r * 0.9
        ax.text(2.0, y, lab, fontsize=7, ha="right", va="center", fontweight="bold")
        for b, (bname, col0) in enumerate(boxes):
            x0 = 2.6 + b * 3.0
            if r == 0:
                ax.text(x0 + 2*cw, 3.95, bname, fontsize=7, ha="center", color="#555", fontweight="bold")
            mot = motif(recs, sp, col0)
            ref = motif(recs, "Homo_sapiens", col0)
            dis = mot != ref
            for j, aa in enumerate(mot):
                x = x0 + j * cw
                ax.add_patch(plt.Rectangle((x - cw/2 + 0.02, y - 0.32), cw - 0.04, 0.64,
                                           facecolor=BG_BAD if dis else BG_OK, edgecolor="none", alpha=0.7))
                ax.text(x, y, aa, fontsize=7, family="monospace", ha="center", va="center",
                        color=COL_NR1 if dis else COL_OK, fontweight="bold")


def main():
    recs = load()
    fig = plt.figure(figsize=(6.5, 8.0))
    gs = gridspec.GridSpec(2, 1, figure=fig, height_ratios=[1.0, 0.32], hspace=0.12,
                           top=0.95, bottom=0.03, left=0.03, right=0.97)
    panel_a(fig.add_subplot(gs[0]), recs)
    panel_b(fig.add_subplot(gs[1]), recs)
    fig.text(0.02, 0.96, "(a)", fontsize=14, fontweight="bold")
    fig.text(0.02, 0.27, "(b)", fontsize=14, fontweight="bold")
    for ext in ("png", "pdf", "svg"):
        fig.savefig(OUT / f"Fig1_gbgc.{ext}", dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(); print("Saved Fig1_gbgc")


if __name__ == "__main__":
    main()
