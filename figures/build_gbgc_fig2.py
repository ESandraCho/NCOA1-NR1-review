#!/usr/bin/env python3
"""
gBGC paper — Figure 2: genomic GC-conversion tract at the NCOA1 NR1 locus.

(a) Per-clade GC at the NR1 region across the 48-species set: ancestrally GC-leaning, with
    reptiles/birds already elevated (~63%) above marsupials (~40%) and placentals (~44%) and
    amphibians (~48%), intensified to an extreme in monotremes (~78%, where alone the motif is
    broken to PGQLP). All values are computed live from the alignment in panel_c.
(b) GC profile across the NCOA1 genomic locus, centered on the NR1 exon, for monotremes
    (platypus + echidna) vs therian (human) and reptile (chicken) references. The NR1 GC peak
    follows an ancestral gradient and is HIGHEST in monotremes: NR1 exon-core GC is 44% in human,
    67% in chicken (reptiles already elevated), 77% echidna, 85% platypus. The therian human is
    the low baseline; the reptile is intermediate.
(c) Fine-grained GC across the NR1 exon: the high-GC tract is a localized peak that roughly
    coincides with the large (~1.3 kb) NR1 coding exon (RefSeq exon 8), strongest in monotremes,
    with the reptile intermediate and the therian low.
(Panel order follows the Results narrative; the panel_* functions are named by content.)

Reads:
  analysis/gbgc_genomic/{platypus,echidna,human,chicken}_NCOA1_genomic.fasta
    (reference loci fetched by fetch_reference_genomic.py, centered on the NR1 exon)
  analysis/outgroups_48sp/NCOA1_48sp_codon_aln.fasta
Outputs: figures/Fig2_gbgc.{png,pdf,svg}
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
from Bio import SeqIO

PAPER = Path(__file__).resolve().parents[1]
OUT = Path(__file__).resolve().parent

# Locate the NR1 exon by motif anchoring (computed from data) rather than hard-coded offsets,
# and orient each locus onto its coding strand. Shares the proven logic with the in-text values.
sys.path.insert(0, str(PAPER / "analysis" / "gbgc_genomic"))
from compute_reported_values import orient_and_find  # noqa: E402

COL_PLAT = "#1565C0"
COL_ECH = "#00897B"
COL_NR1 = "#E53935"
COL_HUM = "#9E9E9E"   # therian reference (grey)
COL_CHK = "#BCAAA4"   # reptile reference (taupe)

# monotremes = the gBGC lineages (foreground); refs = therian + reptile controls
# (name, colour, linestyle): human dashed, chicken densely dotted, to keep the two
# muted controls visually distinct.
MONOTREMES = [("platypus", COL_PLAT, "-"), ("echidna", COL_ECH, "-")]
REFERENCES = [("human", COL_HUM, "--"), ("chicken", COL_CHK, (0, (1, 1)))]


def gc(seq):
    seq = seq.upper()
    return 100 * (seq.count("G") + seq.count("C")) / len(seq) if seq else 0


HALF = 13000  # half-window (bp) around the NR1 exon -> ~26 kb, comparable across species


def load_genomic(name):
    """Return (oriented locus windowed to ~26 kb around NR1, offset of the NR1 motif within it).

    The exon is located by motif anchoring (orient_and_find) — not a hard-coded offset — so the
    profile is always centred on the true NR1 exon even though the fetched loci differ in strand
    and length (e.g. the echidna locus is ~49 kb with NR1 far from its centre)."""
    raw = str(SeqIO.read(PAPER / "analysis" / "gbgc_genomic" / f"{name}_NCOA1_genomic.fasta",
                         "fasta").seq).upper()
    seq, off = orient_and_find(raw)
    if off is None:
        raise RuntimeError(f"NR1 exon not located in {name} genomic locus")
    s0 = max(0, off - HALF)
    return seq[s0:off + HALF], off - s0


def _profile(g, win=300, step=100):
    xs, ys = [], []
    for w in range(0, len(g) - win, step):
        xs.append((w + win / 2) / 1000.0)  # kb
        ys.append(gc(g[w:w + win]))
    ys = np.array(ys)
    sm = np.convolve(ys, np.ones(3) / 3, mode="same")  # light smoothing
    return np.array(xs), sm, ys


def panel_a(ax):
    """GC profile across the locus. The NR1 GC peak follows an ancestral gradient and is
    strongest in monotremes: the reptile (chicken) control is intermediate and the therian
    (human) lowest at NR1, but all show a local peak there."""
    # references first (thin, in the background; human dashed, chicken densely dotted)
    for name, col, ls in REFERENCES:
        g, exon = load_genomic(name)
        # center each profile on its own NR1 exon so positions are comparable
        xs, sm, ys = _profile(g)
        xs = xs - exon / 1000.0
        ax.plot(xs, sm, color=col, lw=0.9, ls=ls, alpha=0.85,
                label=f"{name} (control)", zorder=2)
    # monotremes on top (solid, with NR1 peak marker)
    for name, col, ls in MONOTREMES:
        g, exon = load_genomic(name)
        xs, sm, ys = _profile(g)
        xs = xs - exon / 1000.0
        ax.plot(xs, sm, color=col, lw=1.2, ls=ls, label=name, alpha=0.95, zorder=4)
        nr1_gc = np.mean([y for x, y in zip(xs, ys) if abs(x) < 0.3])
        ax.plot(0, nr1_gc, "o", color=COL_NR1, ms=6, zorder=5)
    ax.axvline(0, color=COL_NR1, ls=":", lw=0.9, alpha=0.7)
    ax.axhline(50, color="#BBBBBB", lw=0.5, ls="--")
    ax.text(0, 96, "NR1 (GC peak\nhighest in monotremes)", color=COL_NR1, fontsize=6.5,
            ha="center", fontweight="bold")
    ax.set_xlabel("Position relative to NR1 exon (kb)", fontsize=8)
    ax.set_ylabel("genomic GC content (%)\n300 bp windows", fontsize=8)
    ax.set_title("Genomic GC profile: NR1 peak strongest in monotremes", fontsize=9, fontweight="bold")
    ax.set_xlim(-12, 12)
    ax.set_ylim(20, 100)
    ax.legend(fontsize=6, loc="lower right", ncol=2)
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(labelsize=7)


# NR1 coding-exon extent (bp relative to the motif start), from RefSeq annotation (exon #8 of
# NCOA1; verified 2026-06-08). The NR1 exon is LARGE (~1320 bp) and the motif sits deep inside it
# — so the GC tract does NOT overrun the exon; it is a localized sub-exonic GC peak.
#   human   NM exon #8: [-799, +522] bp (1321 bp)   platypus exon #8: [-853, +466] bp (1319 bp)
NR1_EXON = {"human": (-0.799, 0.522), "platypus": (-0.853, 0.466)}


def panel_b(ax):
    """Fine GC across the NR1 exon. The high-GC tract is a localized peak that roughly coincides
    with the large (~1.3 kb) NR1 coding exon, strongest in monotremes; the reptile control is
    intermediate and the therian low."""
    win = 150
    offs = list(range(-1500, 1600, 75))
    # references: human dashed line, chicken densely dotted line + small dot markers
    for name, col, ls in REFERENCES:
        g, pos = load_genomic(name)
        xs = [o / 1000.0 for o in offs]
        ys = [gc(g[pos + o: pos + o + win]) for o in offs]
        marker = "." if name == "chicken" else None
        ax.plot(xs, ys, color=col, lw=0.9, ls=ls, marker=marker, markersize=3,
                markevery=1, label=f"{name} (control)", alpha=0.85, zorder=2)
    # monotremes (the tract)
    for name, col, ls in MONOTREMES:
        g, pos = load_genomic(name)
        ys = [gc(g[pos + o: pos + o + win]) for o in offs]
        ax.plot([o / 1000.0 for o in offs], ys, color=col, lw=1.2, ls=ls, label=name, alpha=0.95, zorder=4)
    ax.axvline(0, color=COL_NR1, ls=":", lw=0.8)
    ax.text(0, 92, "NR1 motif", color=COL_NR1, fontsize=7, ha="center", fontweight="bold")
    # mark the NR1 coding exon (RefSeq exon #8, ~1.3 kb) as a bar; the tract is contained within it.
    ex_lo, ex_hi = NR1_EXON["human"]
    y_ex = 25.0
    ax.plot([ex_lo, ex_hi], [y_ex, y_ex], color="#212121", lw=5.0, solid_capstyle="butt", zorder=6)
    for xb in (ex_lo, ex_hi):
        ax.plot([xb, xb], [y_ex - 1.3, y_ex + 1.3], color="#212121", lw=1.0, zorder=6)
    ax.text((ex_lo + ex_hi) / 2, y_ex - 3.0, "NR1 coding exon (~1.3 kb, RefSeq exon 8)",
            fontsize=6, ha="center", va="top", color="#212121")
    ax.text(-1.45, 90, "localized GC peak\nwithin the NR1 exon\n(strongest in monotremes)", fontsize=6.5,
            ha="left", va="top", color="#8D6E00", style="italic")
    ax.set_xlabel("Distance from NR1 motif (kb)", fontsize=8)
    ax.set_ylabel("genomic GC content (%)\n150 bp windows", fontsize=8)
    ax.set_title("Localized GC peak at the NR1 exon (reptile intermediate, extreme in monotremes)",
                 fontsize=8.5, fontweight="bold")
    ax.set_ylim(20, 95)
    ax.legend(fontsize=6, loc="lower right", ncol=2)
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(labelsize=7)


def panel_c(ax):
    """Per-clade NR1-region GC across the 48-species set."""
    recs = {r.id.replace("_NCOA1", ""): str(r.seq)
            for r in SeqIO.parse(PAPER / "analysis" / "outgroups_48sp" / "NCOA1_48sp_codon_aln.fasta", "fasta")}
    # locate NR1 codon via human
    full = {'CTG':'L','CTC':'L','CTA':'L','TTG':'L','TTA':'L','GTG':'V','GTC':'V','GTA':'V','GTT':'V',
            'CAG':'Q','CAA':'Q','CCT':'P','CCC':'P','CCA':'P','CCG':'P','GGG':'G','GGC':'G','GGA':'G','GGT':'G'}
    hum = recs["Homo_sapiens"]
    # translate to find LVQLL codon index
    from Bio.Seq import Seq
    hp = str(Seq(hum.replace("-", "")).translate())
    mi = hp.find("LVQLL")
    # map ungapped motif index to aligned codon col
    ung = 0; col = None
    for c in range(len(hum) // 3):
        cod = hum[c*3:c*3+3]
        if cod != "---":
            if ung == mi: col = c; break
            ung += 1
    lo, hi = col - 15, col + 10

    def region_gc(sp):
        s = recs[sp][lo*3:hi*3].replace("-", "")
        return gc(s)

    # Classify EVERY species in the alignment by clade (explicit membership sets so the
    # assignment is reproducible and reviewer-auditable). Marsupials + placentals shown
    # separately rather than lumped as "Therians" — uses all 48 species, none dropped.
    MONOTREMES = {"Ornithorhynchus_anatinus", "Tachyglossus_aculeatus"}
    MARSUPIALS = {"Monodelphis_domestica", "Sarcophilus_harrisii", "Dromiciops_gliroides",
                  "Petaurus_breviceps_papuanus", "Phascolarctos_cinereus", "Trichosurus_vulpecula",
                  "Vombatus_ursinus"}
    AMPHIBIANS = {"Rana_temporaria", "Bufo_bufo", "Geotrypetes_seraphini", "Xenopus_tropicalis"}
    REPT_BIRDS = {"Chelonia_mydas", "Gallus_gallus", "Anolis_carolinensis", "Alligator_mississippiensis",
                  "Crocodylus_porosus", "Python_bivittatus", "Podarcis_muralis", "Sceloporus_undulatus",
                  "Chrysemys_picta", "Pelodiscus_sinensis", "Thamnophis_elegans"}
    COELACANTH = {"Latimeria_chalumnae"}

    def clade_of(s):
        if s in COELACANTH: return "Coelacanth"
        if s in AMPHIBIANS: return "Amphibians"
        if s in REPT_BIRDS: return "Reptiles/birds"
        if s in MONOTREMES: return "Monotremes"
        if s in MARSUPIALS: return "Marsupials"
        return "Placentals"  # everything else = placental mammal

    # display order: basal -> derived
    order = ["Coelacanth", "Amphibians", "Reptiles/birds", "Monotremes", "Marsupials", "Placentals"]
    cmap = {"Coelacanth": "#90A4AE", "Amphibians": "#8BC34A", "Reptiles/birds": "#7E57C2",
            "Monotremes": COL_NR1, "Marsupials": "#26A69A", "Placentals": "#42A5F5"}
    from collections import defaultdict
    groups = defaultdict(list)
    for sp in recs:
        groups[clade_of(sp)].append(region_gc(sp))

    labels, means, errs, cols = [], [], [], []
    for cl in order:
        vals = groups[cl]
        labels.append(f"{cl}\n(n={len(vals)})")
        means.append(np.mean(vals))
        errs.append(np.std(vals) if len(vals) > 1 else 0)
        cols.append(cmap[cl])
    x = range(len(labels))
    ax.bar(x, means, yerr=errs, color=cols, alpha=0.85, edgecolor="#333", linewidth=0.4,
           error_kw=dict(lw=0.7, capsize=2, capthick=0.7, ecolor="#555"))
    for i, m in enumerate(means):
        ax.text(i, m + errs[i] + 1.5, f"{m:.0f}%", ha="center", fontsize=7, fontweight="bold")
    ax.set_xticks(list(x))
    ax.set_xticklabels(labels, fontsize=6.5)
    ax.set_ylabel("coding GC (%)\nNR1 codon window", fontsize=8)
    ax.set_title("NR1 region GC by clade, all 48 species (motif intact except monotremes)",
                 fontsize=8.5, fontweight="bold")
    ax.axhline(50, color="#BBBBBB", lw=0.5, ls="--")
    ax.set_ylim(0, 95)
    mono_i = order.index("Monotremes")
    ax.annotate("PGQLP\n(motif broken)", xy=(mono_i, means[mono_i] + errs[mono_i]),
                xytext=(mono_i + 0.62, 90), fontsize=6.5, ha="center", color=COL_NR1,
                fontweight="bold",
                arrowprops=dict(arrowstyle="->", color=COL_NR1, lw=0.7))
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(labelsize=7)


def main():
    fig = plt.figure(figsize=(7.0, 8.5))
    gs = gridspec.GridSpec(3, 1, figure=fig, hspace=0.45, top=0.95, bottom=0.07, left=0.12, right=0.95)
    # Panel order follows the Results narrative (cited a→b→c):
    #   a = clade GC gradient (overview), b = genomic locus profile, c = exon-boundary tract.
    # (The panel_* functions are named by content, not position.)
    panel_c(fig.add_subplot(gs[0]))   # a: clade gradient
    panel_a(fig.add_subplot(gs[1]))   # b: locus GC profile
    panel_b(fig.add_subplot(gs[2]))   # c: exon-boundary tract
    fig.text(0.02, 0.96, "(a)", fontsize=14, fontweight="bold")
    fig.text(0.02, 0.64, "(b)", fontsize=14, fontweight="bold")
    fig.text(0.02, 0.32, "(c)", fontsize=14, fontweight="bold")
    for ext in ("png", "pdf", "svg"):
        fig.savefig(OUT / f"Fig2_gbgc.{ext}", dpi=300, bbox_inches="tight", facecolor="white")
    plt.close()
    print("Saved Fig2_gbgc")


if __name__ == "__main__":
    main()
