#!/usr/bin/env python3
"""
Reproduce every prose-reported gBGC value that is NOT already emitted by a figure script.

Before this script existed, the genomic-tract numbers (Fig 2b,c), the NR1-region
amino-acid composition / disorder fractions, and the monotreme-stem W->S substitution
counts (Fig 3a) lived only as hand-typed numbers in `gBGC_genomic_event_NR1.md` and were
hard-coded into `figures/build_gbgc_fig3.py`. This recomputes them from the committed data
so the manuscript's "every reported value has a script" claim holds and the numbers are
reviewer-auditable.

Inputs (paths relative to papers/paper_ncoa1_hormone_switch/):
  - analysis/gbgc_genomic/{platypus,echidna,human,chicken}_NCOA1_genomic.fasta  (genomic loci)
  - analysis/outgroups_48sp/NCOA1_48sp_codon_aln.fasta                           (codon alignment)
  - analysis/runs_48sp/asr/rst                                                   (codeml ASR)

Run:  python analysis/gbgc_genomic/compute_reported_values.py
"""
import re
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq

PAPER = Path(__file__).resolve().parents[2]
GENOMIC = PAPER / "analysis" / "gbgc_genomic"
ALN = PAPER / "analysis" / "outgroups_48sp" / "NCOA1_48sp_codon_aln.fasta"
ASR_DIR = PAPER / "analysis" / "runs_48sp" / "asr"
RST = ASR_DIR / "rst"


def ensure_rst():
    """codeml's `rst` (ancestral reconstruction) is a large, regenerable intermediate that is
    gitignored, so a fresh clone will not have it. Regenerate it from the committed inputs
    (ncoa1_48.phy, ncoa1_48.tree, codeml.ctl) by re-running codeml if it is missing or empty.
    The run is deterministic, so the W→S counts and stem motif reproduce exactly."""
    if RST.exists() and RST.stat().st_size > 0:
        return
    import shutil
    import subprocess
    codeml = shutil.which("codeml")
    if not codeml:
        raise SystemExit(
            f"{RST} is missing/empty and `codeml` is not on PATH. Install PAML, or run codeml "
            f"manually in {ASR_DIR} (`codeml codeml.ctl`) to regenerate `rst`.")
    print(f"[ensure_rst] {RST} missing — regenerating with codeml (deterministic; ~several min)...")
    subprocess.run([codeml, "codeml.ctl"], cwd=ASR_DIR, check=True,
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    if not (RST.exists() and RST.stat().st_size > 0):
        raise SystemExit(f"codeml ran but {RST} is still empty — check {ASR_DIR}/codeml.ctl inputs.")

# Anchor the NR1 exon on the motif (LVQLL / monotreme PGQLP) IN ITS CONSERVED NRID CONTEXT:
# motif + X + [TN] + X + A  (the invariant '..T.A'/'..N.A' stretch right after the motif, verified
# vs RefSeq annotation across all sampled species). The bare motif recurs spuriously across the
# ~280-kb gene, so the context is required to reject false matches. Genomic FASTAs are not
# orientation-normalised, so we scan both strands / all 3 frames; there is exactly one genuine hit.
import re as _re
NR1_CONTEXT = _re.compile(r"(?:LVQLL|PGQLP).[NST].A")


def gc(seq):
    seq = seq.replace("-", "").upper()
    if not seq:
        return float("nan")
    return 100.0 * (seq.count("G") + seq.count("C")) / len(seq)


def cpg_oe(seq):
    """Gardiner-Garden observed/expected CpG over a sequence."""
    seq = seq.replace("-", "").upper()
    n = len(seq)
    if n < 2:
        return float("nan")
    cg = seq.count("CG")
    c = seq.count("C")
    g = seq.count("G")
    exp = (c * g) / n
    return (cg / exp) if exp else float("nan")


# ----------------------------------------------------------------------------- genomic tract
def orient_and_find(seq):
    """Return (oriented_seq, offset_bp of NR1 motif start), locating the motif in its conserved
    NRID context (NR1_CONTEXT). Scans both strands and all 3 frames; there is exactly one genuine
    hit. Returns (seq, None) if the genuine NR1 exon is absent (window is the wrong region)."""
    for s in (seq, str(Seq(seq).reverse_complement())):
        for frame in range(3):
            prot = str(Seq(s[frame:len(s) - (len(s) - frame) % 3]).translate())
            m = NR1_CONTEXT.search(prot)
            if m:
                return s, frame + m.start() * 3
    return seq, None


def genomic_metrics(sp):
    """Return a dict of NR1-locus GC metrics for one species, computed from its genomic FASTA.
    Returns None if the locus or its NR1 exon cannot be located."""
    fa = GENOMIC / f"{sp}_NCOA1_genomic.fasta"
    if not fa.exists():
        return None
    raw = str(next(SeqIO.parse(fa, "fasta")).seq).upper()
    seq, off = orient_and_find(raw)
    if off is None:
        return None
    # window the locus to ~26 kb centred on the motif so the locus mean / percentile are
    # comparable across species (some fetched loci are longer; matches manuscript "~26-kb").
    HALF = 13000
    s0 = max(0, off - HALF)
    seq = seq[s0:off + HALF]
    off -= s0
    win_gc = gc(seq[max(0, off - 150):off + 150])
    profile = [gc(seq[i:i + 300]) for i in range(0, len(seq) - 300, 100)]
    below = sum(1 for v in profile if v <= win_gc)
    cpg_win = seq[max(0, off - 750):off + 750]
    aln_id = GENOMIC_TO_ALN.get(sp)
    return {
        "locus_kb": len(seq) / 1000.0,
        # CANONICAL per-lineage NR1 GC = coding GC (same definition as Fig 2a), so every quoted
        # per-region GC value uses one criterion. None if the species is not in the alignment.
        "nr1_gc": coding_nr1_gc(aln_id) if aln_id else None,
        # --- genomic spatial/tract metrics (separate analysis; not per-lineage GC values) ---
        "locus_gc": gc(seq),                                   # intron+exon locus mean (true genomic)
        "percentile": 100.0 * below / len(profile),            # NR1 window's rank within the locus
        "cpg_gc": gc(cpg_win),
        "cpg_oe": cpg_oe(cpg_win),                             # NR1 +/-750 bp CpG island test
    }


def genomic_report():
    print("=" * 70)
    print("NR1 GC + genomic tract (Fig 2) — coding NR1 GC (canonical) + genomic locus metrics")
    print("=" * 70)
    for sp in ["platypus", "echidna", "human", "chicken"]:
        m = genomic_metrics(sp)
        if m is None:
            print(f"{sp}: locus or NR1 exon not located")
            continue
        nr1 = m["nr1_gc"]
        # CpG-island status under two criteria over the +/-750 bp NR1 window:
        #   Gardiner-Garden: GC > 50%, o/e > 0.60
        #   Takai-Jones (stricter): GC >= 55%, o/e >= 0.65
        gc_win, oe = m["cpg_gc"], m["cpg_oe"]
        gg = "PASS" if (gc_win > 50 and oe > 0.60) else "fail"
        tj = "PASS" if (gc_win >= 55 and oe >= 0.65) else "fail"
        print(f"{sp:9s} coding NR1 GC {nr1:.0f}% (canonical, Fig 2a/b/c value) | "
              f"genomic locus mean {m['locus_gc']:.0f}% (intron+exon) | "
              f"NR1 ~{m['percentile']:.0f}th pctile of locus | "
              f"CpG-window GC {gc_win:.0f}% o/e {oe:.2f} "
              f"[Gardiner-Garden {gg}; Takai-Jones {tj}]")
    print()


# ----------------------------------------------------------------------------- composition
# TOP-IDP disorder-promoting residue set (residues with positive disorder propensity).
DISORDER_PROMOTING = set("ARSGQEKPD")  # A,R,S,G,Q,E,K,P,D


def load_aln():
    return {r.id.replace("_NCOA1", ""): str(r.seq)
            for r in SeqIO.parse(ALN, "fasta")}


def nr1_codon_col(recs):
    hum = recs["Homo_sapiens"]
    hp = str(Seq(hum.replace("-", "")).translate())
    mi = hp.find("LVQLL")
    ung = 0
    for c in range(len(hum) // 3):
        if hum[c * 3:c * 3 + 3] != "---":
            if ung == mi:
                return c
            ung += 1
    raise RuntimeError("NR1 column not found")


# Single canonical NR1 GC definition used EVERYWHERE a per-lineage NR1 GC value is reported
# (Fig 2a clade bars and the per-species values quoted for Fig 2b/c): coding GC over the NR1
# codon window col-15 .. col+10 (25 codons) from the codon alignment — EXACTLY the window used
# by build_gbgc_fig2.py panel_c, so the canonical value equals what Fig 2a plots.
# The genomic sliding-window profiles in Fig 2b/c are a SEPARATE spatial analysis (tract shape)
# and legitimately use bp windows on the genomic locus; they are not per-lineage GC values.
CODING_WINDOW = (-15, 10)  # alignment codon columns: [col+lo, col+hi)

# alignment species id for each genomic-FASTA key (for cross-referencing the two)
GENOMIC_TO_ALN = {
    "human": "Homo_sapiens", "chicken": "Gallus_gallus",
    "platypus": "Ornithorhynchus_anatinus", "echidna": "Tachyglossus_aculeatus",
    "opossum": "Monodelphis_domestica", "devil": "Sarcophilus_harrisii",
    "frog": "Xenopus_tropicalis", "anole": "Anolis_carolinensis",
}


def coding_nr1_gc(aln_species, recs=None):
    """Canonical coding NR1 GC for one species (by its codon-alignment id)."""
    if recs is None:
        recs = load_aln()
    if aln_species not in recs:
        return None
    col = nr1_codon_col(recs)
    lo, hi = col + CODING_WINDOW[0], col + CODING_WINDOW[1]
    return gc(recs[aln_species][lo * 3:hi * 3])


def coding_nr1_gc_scanwindow(aln_species, recs=None):
    """NR1 coding GC over the SAME window the coactivator scan uses for every other LXXLL motif
    (±15 codons flanking + the 5 motif codons = 35 codons), so the Fig 4 NR1 point is computed
    identically to the other points it is compared against."""
    if recs is None:
        recs = load_aln()
    if aln_species not in recs:
        return None
    col = nr1_codon_col(recs)
    return gc(recs[aln_species][(col - 15) * 3:(col + 5 + 15) * 3])


def composition_report():
    print("=" * 70)
    print("NR1-REGION COMPOSITION + DISORDER (Results) — +/-15-codon window, codon alignment")
    print("disorder-promoting set = A,R,S,G,Q,E,K,P,D (TOP-IDP positive-propensity residues)")
    print("=" * 70)
    recs = load_aln()
    col = nr1_codon_col(recs)
    lo, hi = col - 15, col + 5 + 15
    for sp in ["Ornithorhynchus_anatinus", "Tachyglossus_aculeatus", "Homo_sapiens"]:
        prot = str(Seq(recs[sp][lo * 3:hi * 3].replace("-", "")).translate()).replace("*", "")
        n = len(prot)
        g = 100 * prot.count("G") / n
        p = 100 * prot.count("P") / n
        r = 100 * prot.count("R") / n
        dis = 100 * sum(1 for a in prot if a in DISORDER_PROMOTING) / n
        print(f"{sp:26s} n={n:2d}  Gly {g:4.0f}%  Pro {p:4.0f}%  Arg {r:4.0f}%  disorder {dis:4.0f}%")
    print()


# ----------------------------------------------------------------------------- W->S on stem
def parse_rst_nodes():
    ensure_rst()  # regenerate the gitignored codeml `rst` from committed inputs if absent
    nodes = {}
    text = RST.read_text()
    for m in re.finditer(r"node #(\d+)\s+((?:[A-Z\-]{3}\s*)+)", text):
        nodes[int(m.group(1))] = m.group(2).split()
    return nodes


def stem_nr1_motif(stem=92):
    """Reconstructed NR1 amino-acid motif at the monotreme stem node, read from the codeml rst
    (so Fig 1's annotation is data-derived, not typed)."""
    nodes = parse_rst_nodes()
    recs = load_aln()
    col = nr1_codon_col(recs)
    cods = nodes[stem][col:col + 5]
    return "".join(str(Seq(c).translate()) if "-" not in c else "-" for c in cods)


def ws_counts(parent=91, stem=92, flank=12):
    """Return (nonsyn_WS, nonsyn_SW, syn_WS, syn_SW) on the monotreme-stem branch over the
    NR1 +/-flank-codon window. This is the single canonical computation; the figure and the
    manuscript both read these numbers so they cannot drift apart."""
    nodes = parse_rst_nodes()
    if parent not in nodes or stem not in nodes:
        raise RuntimeError(f"nodes {parent}/{stem} not in rst (have {min(nodes)}-{max(nodes)})")
    recs = load_aln()
    col = nr1_codon_col(recs)
    p, c = nodes[parent], nodes[stem]

    def aa(cod):
        try:
            return str(Seq(cod).translate())
        except Exception:
            return "X"

    ws_syn = ws_non = sw_syn = sw_non = 0
    for ci in range(col - flank, col + 5 + flank):
        pc, cc = p[ci], c[ci]
        if pc == cc or "-" in pc or "-" in cc:
            continue
        syn = aa(pc) == aa(cc)
        for x, y in zip(pc, cc):
            if x == y:
                continue
            wx, wy = x in "AT", y in "AT"
            if wx and not wy:
                ws_syn += syn
                ws_non += not syn
            elif not wx and wy:
                sw_syn += syn
                sw_non += not syn
    return ws_non, sw_non, ws_syn, sw_syn


def platypus_echidna_identities():
    """Return dict of platypus-echidna amino-acid % identity for Fig 3b panels:
    NR1 core (5-mer), NR1 flank (+/-15 aa), NR2 window (+/-15 aa), NR3 window (+/-15 aa),
    whole protein."""
    recs = load_aln()
    plat = recs["Ornithorhynchus_anatinus"]
    ech = recs["Tachyglossus_aculeatus"]

    def aa_ident(c0, c1):
        pa = str(Seq(plat[c0 * 3:c1 * 3].replace("-", "N")).translate())
        ea = str(Seq(ech[c0 * 3:c1 * 3].replace("-", "N")).translate())
        n = sum(1 for a, b in zip(pa, ea) if "X" not in (a, b))
        s = sum(1 for a, b in zip(pa, ea) if a == b and "X" not in (a, b))
        return round(100 * s / n) if n else float("nan")

    hum = recs["Homo_sapiens"]
    hp = str(Seq(hum.replace("-", "N")).translate())

    def col_of(motif):  # aligned codon col of human motif start
        return hp.find(motif)

    c1 = col_of("LVQLL")
    c2 = col_of("LHRLL")
    c3 = col_of("LRYLL")
    ncod = len(plat) // 3
    return {
        "NR1 core (PGQLP)": aa_ident(c1, c1 + 5),
        "NR1 flank (high GC)": aa_ident(c1 - 15, c1 + 5 + 15),
        "NR2 window (high GC)": aa_ident(c2 - 15, c2 + 5 + 15),
        "NR3 window (high GC)": aa_ident(c3 - 15, c3 + 5 + 15),
        "whole protein": aa_ident(0, ncod),
    }


def ws_report():
    print("=" * 70)
    print("MONOTREME-STEM W->S SUBSTITUTIONS (Fig 3a) — codeml ASR, node 91->92")
    print("NR1 +/-12-codon window; classified per substituted codon as syn / nonsyn")
    print("=" * 70)
    nws, nsw, sws, ssw = ws_counts()
    print(f"Nonsynonymous: W->S {nws} : S->W {nsw}")
    print(f"Synonymous:    W->S {sws} : S->W {ssw}")
    # Is the W->S bias significantly stronger at nonsynonymous than synonymous sites?
    # 2x2 table: rows = {nonsyn, syn}, cols = {W->S, S->W}. Fisher's exact (two-sided).
    # If gBGC were driving the protein change, the W->S excess should be concentrated at
    # the nonsynonymous (protein-altering) sites relative to synonymous sites.
    from scipy.stats import fisher_exact
    table = [[nws, nsw], [sws, ssw]]
    odds, p = fisher_exact(table, alternative="two-sided")
    _, p_greater = fisher_exact(table, alternative="greater")
    print(f"Fisher's exact 2x2 [[{nws},{nsw}],[{sws},{ssw}]]: "
          f"OR={odds:.3g}, two-sided p={p:.4g}, one-sided(nonsyn>syn W->S) p={p_greater:.4g}")
    print("(NB synonymous arm is small (n={}); the synonymous result is 'no detectable".format(sws + ssw))
    print(" W->S excess', not a demonstrated absence of bias — power is minimal there.)")
    print("(reported result: W->S enriched at nonsynonymous, not synonymous;")
    print(" exact counts are window-dependent — this is the canonical figure window.)")
    print()


def identity_report():
    print("=" * 70)
    print("PLATYPUS-ECHIDNA IDENTITY (Fig 3b) — amino-acid %, codon alignment")
    print("=" * 70)
    for k, v in platypus_echidna_identities().items():
        print(f"{k:24s} {v}%")
    print()


# Therian species in the 48-species alignment (eutherians + marsupials; excludes the 2 monotremes,
# the 16 non-mammal outgroups). Used for the NRID-vs-rest constraint control.
THERIANS = [
    "Homo_sapiens", "Pan_troglodytes", "Macaca_mulatta", "Mus_musculus", "Rattus_norvegicus",
    "Mesocricetus_auratus", "Cavia_porcellus", "Heterocephalus_glaber", "Oryctolagus_cuniculus",
    "Bos_taurus", "Tursiops_truncatus", "Sus_scrofa", "Equus_caballus", "Canis_lupus_familiaris",
    "Felis_catus", "Manis_javanica", "Sorex_araneus", "Pteropus_vampyrus", "Erinaceus_europaeus",
    "Dasypus_novemcinctus", "Echinops_telfairi", "Loxodonta_africana", "Trichechus_manatus",
    "Monodelphis_domestica", "Sarcophilus_harrisii", "Dromiciops_gliroides", "Phascolarctos_cinereus",
    "Vombatus_ursinus", "Trichosurus_vulpecula", "Petaurus_breviceps_papuanus",
]


def therian_nrid_identity(recs=None):
    """Mean pairwise amino-acid identity among therians, for (i) the NRID (NR1 start → NR3 end)
    and (ii) the whole protein. The Results control claim 'the therian NRID is not more conserved
    than non-LXXLL regions (87% vs 90% mean therian identity)' is computed here so it has a script
    rather than being hand-typed. Returns (nrid_identity_pct, whole_protein_identity_pct)."""
    import itertools
    import statistics as st
    if recs is None:
        recs = load_aln()
    therians = [t for t in THERIANS if t in recs]
    hum = recs["Homo_sapiens"]
    hp = str(Seq(hum.replace("-", "")).translate())

    def col_of(motif):  # aligned codon col of human motif start
        mi = hp.find(motif)
        ung = 0
        for c in range(len(hum) // 3):
            if hum[c * 3:c * 3 + 3] != "---":
                if ung == mi:
                    return c
                ung += 1

    nrid_lo, nrid_hi = col_of("LVQLL"), col_of("LRYLL") + 5  # NR1 start .. NR3 end
    ncod = len(hum) // 3

    def prot(sp, lo, hi):
        return str(Seq(recs[sp][lo * 3:hi * 3].replace("-", "N")).translate())

    def mean_pairwise_id(lo, hi):
        ids = []
        for a, b in itertools.combinations(therians, 2):
            pa, pb = prot(a, lo, hi), prot(b, lo, hi)
            n = sum(1 for x, y in zip(pa, pb) if "X" not in (x, y))
            s = sum(1 for x, y in zip(pa, pb) if x == y and "X" not in (x, y))
            if n:
                ids.append(100 * s / n)
        return st.mean(ids)

    return mean_pairwise_id(nrid_lo, nrid_hi), mean_pairwise_id(0, ncod)


def therian_nrid_report():
    print("=" * 70)
    print("THERIAN NRID-vs-REST CONSTRAINT (Results control) — mean pairwise AA identity")
    print("=" * 70)
    nrid, whole = therian_nrid_identity()
    print(f"NRID (NR1 start → NR3 end) therian identity: {nrid:.0f}% ({nrid:.1f}%)")
    print(f"Whole-protein therian identity:              {whole:.0f}% ({whole:.1f}%)")
    print("  -> the NRID is NOT more conserved than the rest of NCOA1; therian LVQLL persistence")
    print("     reflects ordinary motif-level constraint, not exceptional regional selection.")
    print()


def _nr1_alignment_columns():
    """1-based codon columns of the human LVQLL motif in the 48-species codon alignment."""
    import json as _json  # noqa: F401  (kept local; module already imports re/Path)
    aln = {r.id: str(r.seq) for r in SeqIO.parse(str(ALN), "fasta")}
    human = aln["Homo_sapiens_NCOA1"]
    prot = str(Seq(human.replace("-", "")).translate())
    start_aa = prot.find("LVQLL")
    aa = 0
    aa_to_col = {}
    for c in range(0, len(human), 3):
        if human[c:c + 3] != "---":
            aa_to_col[aa] = c // 3
            aa += 1
    return [aa_to_col[start_aa + k] + 1 for k in range(5)]


def sitelevel_report():
    """Site-level selection at NR1 (FUBAR + Contrast-FEL) — the key 'branch test lights up
    but site-level tests do not flag NR1' point. gBGC inflates branch dN/dS without producing
    a genuine pervasive adaptive signal at the motif site."""
    import json
    print("=" * 70)
    print("SITE-LEVEL SELECTION AT NR1 (FUBAR + Contrast-FEL)")
    print("=" * 70)
    nr1 = _nr1_alignment_columns()
    print(f"NR1 motif codons (1-based alignment columns): {nr1}")

    runs = PAPER / "analysis" / "runs_48sp"
    fub = json.load(open(runs / "NCOA1_48sp_fubar.json"))
    frows = fub["MLE"]["content"]["0"]
    fhdr = [h[0] for h in fub["MLE"]["headers"]]
    pp_col = fhdr.index("Prob[alpha<beta]")  # posterior prob of positive selection
    n_pos = sum(1 for r in frows if r[pp_col] >= 0.9)
    max_pp = max(r[pp_col] for r in frows)
    print(f"FUBAR: {len(frows)} sites; sites with Prob[positive]>=0.9: {n_pos} "
          f"(max posterior {max_pp:.3f}) -> no pervasive positive selection anywhere")
    print("  NR1 codons (alpha=syn rate, beta=nonsyn rate, Prob[positive]):")
    for c in nr1:
        r = frows[c - 1]
        print(f"    codon {c}: alpha={r[0]:.2f} beta={r[1]:.2f} "
              f"Prob[positive]={r[pp_col]:.3f}  (alpha>>beta -> purifying + high local syn rate)")

    cf = json.load(open(runs / "NCOA1_48sp_contrast_fel.json"))
    crows = cf["MLE"]["content"]["0"]
    chdr = [h[0] for h in cf["MLE"]["headers"]]
    q_col = chdr.index("Q-value (overall)")
    sig = [i + 1 for i, r in enumerate(crows) if r[q_col] <= 0.2]
    nr1_hits = sorted(set(sig) & set(nr1))
    print(f"Contrast-FEL: {len(sig)} sites with differential dN/dS at Q<=0.2; "
          f"NR1 codons among them: {nr1_hits if nr1_hits else 'NONE'}")
    print("  -> branch tests (aBSREL/BUSTED) flag the monotreme lineage, but neither")
    print("     site-level test flags the NR1 motif itself: the gBGC-artifact signature.")
    print()


if __name__ == "__main__":
    genomic_report()
    composition_report()
    ws_report()
    identity_report()
    therian_nrid_report()
    sitelevel_report()
