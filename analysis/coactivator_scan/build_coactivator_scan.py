#!/usr/bin/env python3
"""
Coactivator LXXLL scan — regenerates coactivator_lxxll_gc.csv for Fig 4.

For every LXXLL ("NR box") motif in nine nuclear-receptor coactivators, compute the local GC
content of the surrounding coding sequence in monotremes, and record whether the motif is intact
or disrupted in monotremes relative to human.

METHOD (documented; this is what Methods / Fig 4 legend describe):
  1. For each gene, the human, platypus and echidna proteins are aligned with MAFFT (L-INS-i),
     the same aligner used throughout the pipeline. (Human↔monotreme coding sequences are not
     index-aligned because of indels, so per-gene alignment is required to identify orthologous
     motifs; a naive same-index lookup mislocates motifs and is not used.)
  2. Each human LXXLL motif (regex L..LL) is located by its alignment columns; the orthologous
     monotreme residues are read from the SAME alignment columns and mapped back to each
     monotreme's ungapped protein/codon index.
  3. Local GC = the GC content of a window of ±15 codons flanking the motif WITH the 5 motif
     codons included (15 + 5 + 15 = 35 codons = 105 nt), measured around the orthologous codons
     in platypus and in echidna. The reported monotreme_local_GC is the MEAN of the two; the
     per-species platypus_GC and echidna_GC are also written (Fig 4 plots the mean with a whisker
     spanning the two — the two agree within ~3% for all motifs except NR1).
  4. Status: a motif is "disrupted" if the orthologous monotreme residues no longer match the
     LXXLL consensus (L at positions 1, 4, 5) in EITHER monotreme; otherwise "intact". (The EP300
     LQNLL case lies in a natively low-complexity Pro/Ala-rich region of p300 and is an alignment
     artifact; it is re-set to intact downstream in build_gbgc_fig4.py, not here.)

NB: NCOA1 NR1 (LVQLL→PGQLP) — the heavily restructured monotreme NR1 region does not align cleanly,
so its scan value is unreliable; build_gbgc_fig4.py adds NR1 with its genomic GC (computed,
mean ~79%) instead.

Requires: mafft on PATH.
Reads:    coactivator_cds.fasta   (human + platypus + echidna CDS for the 9 genes)
Writes:   coactivator_lxxll_gc.csv
"""
import csv
import re
import subprocess
import tempfile
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq

HERE = Path(__file__).resolve().parent
GENES = ["NCOA1", "NCOA2", "NCOA3", "NCOA6", "PPARGC1A", "MED1", "NRIP1", "EP300", "NCOA4"]
FLANK_CODONS = 15  # ±15 codons around the motif, motif included


def gc(s):
    s = s.upper()
    return 100 * (s.count("G") + s.count("C")) / len(s) if s else 0.0


def translate(cds):
    return str(Seq(cds[:len(cds) // 3 * 3]).translate())


def local_gc(cds, codon_index):
    """GC% over ±FLANK_CODONS codons around the motif (motif included)."""
    lo = max(0, (codon_index - FLANK_CODONS) * 3)
    hi = min(len(cds), (codon_index + 5 + FLANK_CODONS) * 3)
    return gc(cds[lo:hi])


def find_lxxll(prot):
    """All LXXLL motifs (L..LL) in a protein, as (0-based aa index, motif)."""
    return [(m.start(), m.group()) for m in re.finditer(r"L..LL", prot)]


def is_lxxll(m):
    return len(m) == 5 and m[0] == "L" and m[3] == "L" and m[4] == "L"


def mafft_align(seqs_by_name):
    """Align proteins with MAFFT L-INS-i. Returns {name: aligned_seq}."""
    with tempfile.NamedTemporaryFile("w", suffix=".fasta", delete=False) as f:
        for name, s in seqs_by_name.items():
            f.write(f">{name}\n{s}\n")
        inp = f.name
    out = subprocess.run(["mafft", "--localpair", "--maxiterate", "1000", "--quiet", inp],
                         capture_output=True, text=True, check=True).stdout
    aln = {}
    name = None
    for line in out.splitlines():
        if line.startswith(">"):
            name = line[1:].strip(); aln[name] = ""
        elif name:
            aln[name] += line.strip()
    Path(inp).unlink(missing_ok=True)
    return aln


def col_to_ungapped(aligned, col):
    """Map an alignment column to the ungapped index in that sequence (or None if a gap)."""
    if aligned[col] == "-":
        return None
    return len(aligned[:col].replace("-", ""))


def ungapped_to_col(aligned, idx):
    """Map an ungapped residue index to its alignment column."""
    seen = 0
    for c, ch in enumerate(aligned):
        if ch != "-":
            if seen == idx:
                return c
            seen += 1
    return None


def main():
    seqs = {r.id: str(r.seq).upper() for r in SeqIO.parse(HERE / "coactivator_cds.fasta", "fasta")}
    out_rows = []
    for gene in GENES:
        hum, plat, ech = (seqs.get(f"{gene}_{sp}") for sp in ("human", "platypus", "echidna"))
        if not (hum and plat and ech):
            continue
        hp, pp, ep = translate(hum), translate(plat), translate(ech)
        aln = mafft_align({"human": hp, "platypus": pp, "echidna": ep})
        ah, ap, ae = aln["human"], aln["platypus"], aln["echidna"]

        for hidx, hmot in find_lxxll(hp):
            col = ungapped_to_col(ah, hidx)          # human motif start column
            if col is None:
                continue
            # orthologous start index in each monotreme (same alignment column)
            pj = col_to_ungapped(ap, col)
            ej = col_to_ungapped(ae, col)
            pmot = pp[pj:pj + 5] if pj is not None else ""
            emot = ep[ej:ej + 5] if ej is not None else ""
            vals = []
            if pj is not None:
                vals.append(local_gc(plat, pj))
            if ej is not None:
                vals.append(local_gc(ech, ej))
            if not vals:
                continue
            mono_gc = round(sum(vals) / len(vals), 1)
            pgc = round(local_gc(plat, pj), 1) if pj is not None else ""
            egc = round(local_gc(ech, ej), 1) if ej is not None else ""
            disrupted = (pmot and not is_lxxll(pmot)) or (emot and not is_lxxll(emot))
            out_rows.append({
                "gene": gene,
                "human_motif": hmot,
                "pos": hidx + 1,
                "monotreme_local_GC": mono_gc,
                "platypus_GC": pgc,
                "echidna_GC": egc,
                "status": "disrupted" if disrupted else "intact",
                "monotreme_motifs": ";".join([m for m in (pmot, emot) if m]),
            })

    with open(HERE / "coactivator_lxxll_gc.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["gene", "human_motif", "pos", "monotreme_local_GC",
                                          "platypus_GC", "echidna_GC", "status", "monotreme_motifs"])
        w.writeheader()
        w.writerows(out_rows)
    print(f"Wrote {len(out_rows)} motif rows (MAFFT-aligned, ±{FLANK_CODONS} codons motif-included, "
          f"platypus/echidna mean)")


if __name__ == "__main__":
    main()
