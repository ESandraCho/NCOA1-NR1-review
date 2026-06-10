#!/usr/bin/env python3
"""
Validate that every *_NCOA1_genomic.fasta window actually contains the genuine NR1 exon.

WHY THIS EXISTS: the NR1 exon is located by motif anchoring. A naive 'SHKL' or bare 'LVQLL'
anchor can match a SPURIOUS site elsewhere in a mis-fetched window, silently centring the GC
profile on the wrong locus. This checker requires the motif IN ITS CONSERVED NRID CONTEXT
(LVQLL / monotreme PGQLP immediately followed by the invariant '..TAEQ' / 'NTAEQ' NRID flank),
which a spurious LVQLL elsewhere in the genome will not satisfy.

A window FAILS if the genuine NR1 exon is absent (window is the wrong genomic region) — in which
case any GC/CpG value computed from it (Fig 2b/2c, Fig S1) is meaningless and must not be used.

Run:   python analysis/gbgc_genomic/validate_genomic_fastas.py
Exit:  0 if all present, 1 if any window is missing the genuine NR1 exon.
"""
import glob
import os
import re
import sys
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq

GENOMIC = Path(__file__).resolve().parent

# The NR1 motif (LVQLL, or monotreme PGQLP) in its conserved NRID downstream signature:
# motif + X + [TN] + X + A  ('..T.A' / '..N.A' — invariant across all sampled species; verified
# vs RefSeq annotation). A spurious LVQLL elsewhere in the genome will not carry this context.
NR1_CONTEXT = re.compile(r"(?:LVQLL|PGQLP).[NST].A")


def has_genuine_nr1(seq):
    """True if a strand/frame translation contains the NR1 motif in conserved NRID context."""
    rc = str(Seq(seq).reverse_complement())
    for s in (seq, rc):
        for f in range(3):
            prot = str(Seq(s[f:len(s) - (len(s) - f) % 3]).translate())
            m = NR1_CONTEXT.search(prot)
            if m:
                return True, prot[m.start():m.start() + 12]
    return False, None


def main():
    failures = []
    for fa in sorted(glob.glob(str(GENOMIC / "*_NCOA1_genomic.fasta"))):
        sp = os.path.basename(fa).replace("_NCOA1_genomic.fasta", "")
        seq = str(next(SeqIO.parse(fa, "fasta")).seq).upper()
        ok, ctx = has_genuine_nr1(seq)
        if ok:
            print(f"  OK      {sp:12s} ({len(seq)} bp)  {ctx}")
        else:
            print(f"  FAIL    {sp:12s} ({len(seq)} bp)  genuine NR1 exon ABSENT — wrong region")
            failures.append(sp)
    if failures:
        print(f"\n{len(failures)} window(s) do NOT contain the genuine NR1 exon: {', '.join(failures)}")
        print("Re-fetch these with analysis/gbgc_genomic/fetch_reference_genomic.py before using "
              "any GC/CpG value derived from them.")
        return 1
    print("\nAll genomic windows contain the genuine NR1 exon.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
