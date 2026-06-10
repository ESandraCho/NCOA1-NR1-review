#!/usr/bin/env python3
"""
Fetch homologous NCOA1 *genomic* loci (intron-containing) for a therian (human) and
reptiles (chicken, anole) so Fig 2 panels a/b can show the monotreme GC peak /
tract-beyond-exon is monotreme-SPECIFIC, not an amniote-wide feature of the locus.

NCOA1 is a large gene (280-340 kb), so we don't fetch the whole locus. We locate the
NR1 exon by its coding motif (the conserved LVQLL...; in monotremes PGQLP) and fetch a
~26 kb window centered on it, matching the monotreme genomic windows. GC profiles are
then directly comparable across clades at the same locus position.

Strategy per species:
  1. esummary(gene) -> chromosome accession + gene span (GenomicInfo).
  2. efetch the full gene region (padded), 6-frame scan for the NR1 peptide to find the
     exon, then re-window +/-13 kb around it.
Outputs: analysis/gbgc_genomic/{human,chicken,anole}_NCOA1_genomic.fasta
"""
import os
import time
from pathlib import Path
from Bio import Entrez, SeqIO
from Bio.Seq import Seq

# NCBI Entrez requires a contact email. Set ENTREZ_EMAIL in your environment before running.
Entrez.email = os.environ.get("ENTREZ_EMAIL", "your.email@example.com")
OUT = Path(__file__).resolve().parent

SPECIES = {
    "human":   "Homo sapiens",
    "chicken": "Gallus gallus",
    "anole":   "Anolis carolinensis",
    # deeper outgroups added 2026-06-08 to anchor the ANCESTRAL CpG-island state
    "turtle":     "Chelonia mydas",
    "crocodile":  "Crocodylus porosus",
    "frog":       "Xenopus tropicalis",
    # marsupial: tests whether CpG-island erosion is therian-wide (PRDM9) vs placental-only
    "opossum":    "Monodelphis domestica",
    "devil":      "Sarcophilus harrisii",
}
HALF = 13000           # half-window (bp) around the NR1 exon -> ~26 kb, matches monotremes
# NR1 peptide context (verified vs RefSeq annotation across all sampled species):
#   human/turtle/devil/opossum/frog  LVQLL-ATTA-EQQ      crocodile/chicken  LVQLL-ASTA-EQQ
#   anole                            LVQLL-ATTA-QEQ      platypus           PGQLP-ANAA-ERR
#                                                        echidna            PGQLP-ANTA-EQR
# The MOTIF ALONE (LVQLL/PGQLP) recurs spuriously across a 280-kb gene, so we REQUIRE the
# conserved NRID downstream signature: motif + X + [NST] + X + A  (position +8 after the motif is
# an invariant 'A' across all 10 sampled species; +6 is N/S/T). A bare SHKL/LVQLL match without it
# is rejected. We scan BOTH strands and all frames and take the genuine match (exactly one).
import re as _re
NR1_CONTEXT = _re.compile(r"(?:LVQLL|PGQLP).[NST].A")


def gene_locus(species):
    h = Entrez.esearch(db="gene", term=f'NCOA1[sym] AND "{species}"[orgn]')
    ids = Entrez.read(h)["IdList"]; h.close()
    if not ids:
        raise RuntimeError("no gene id")
    h = Entrez.esummary(db="gene", id=ids[0])
    doc = Entrez.read(h, validate=False)["DocumentSummarySet"]["DocumentSummary"][0]; h.close()
    g = doc["GenomicInfo"][0]
    acc = g["ChrAccVer"]
    a, b = int(g["ChrStart"]), int(g["ChrStop"])
    start, stop = min(a, b), max(a, b)
    strand = "minus" if a > b else "plus"
    return acc, start, stop, strand


def efetch_region(acc, s, e, strand):
    strand_n = 2 if strand == "minus" else 1
    h = Entrez.efetch(db="nuccore", id=acc, rettype="fasta", retmode="text",
                      seq_start=s + 1, seq_stop=e + 1, strand=strand_n)
    rec = SeqIO.read(h, "fasta"); h.close()
    return str(rec.seq).upper()


def find_nr1(seq):
    """Return (oriented_seq, offset_bp) of the genuine NR1 exon — the motif in its conserved NRID
    context (NR1_CONTEXT). Scans both strands and all frames; a bare motif match without the NRID
    context is ignored. Returns (seq, None) if not found."""
    for s in (seq, str(Seq(seq).reverse_complement())):
        for frame in range(3):
            prot = str(Seq(s[frame:len(s) - (len(s) - frame) % 3]).translate())
            m = NR1_CONTEXT.search(prot)
            if m:
                return s, frame + m.start() * 3
    return seq, None


def main():
    for key, sci in SPECIES.items():
        try:
            acc, gs, ge, strand = gene_locus(sci)
            print(f"{key:8s} {sci:22s} {acc}:{gs}-{ge} ({strand}) span={ge-gs}")
            # fetch the whole gene on its coding strand, scan BOTH strands for the genuine NR1
            full = efetch_region(acc, gs, ge, strand)
            full, off = find_nr1(full)
            if off is None:
                print(f"         NR1 NOT FOUND in {key} (genuine NRID context absent) — SKIPPING")
                continue
            print(f"         NR1 exon at offset {off} bp within gene region")
            lo = max(0, off - HALF)
            hi = min(len(full), off + HALF)
            sub = full[lo:hi]
            # the NR1 exon should now sit at off-lo within `sub`
            rec_off = off - lo
            rec = SeqIO.read(  # build a SeqRecord cheaply
                __import__("io").StringIO(f">{key}_NCOA1_locus {acc} window NR1@{rec_off}\n{sub}\n"),
                "fasta")
            SeqIO.write(rec, OUT / f"{key}_NCOA1_genomic.fasta", "fasta")
            print(f"         -> wrote {key}_NCOA1_genomic.fasta ({len(sub)} bp), NR1 at {rec_off}")
        except Exception as ex:
            print(f"{key}: FAILED — {ex}")
        time.sleep(0.5)


if __name__ == "__main__":
    main()
