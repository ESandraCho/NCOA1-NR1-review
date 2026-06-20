# NCOA1 NR1 — data and reproducibility package

Data and analysis code to reproduce every figure and reported value for the study of the
NCOA1 NR-box 1 (NR1) substitution in monotremes and its origin by GC-biased gene conversion.

## Layout

```
analysis/
  outgroups_48sp/    48-species NCOA1 codon + protein alignments and ML tree
  runs_48sp/         HyPhy selection results (aBSREL, BUSTED, Contrast-FEL, FUBAR)
                     and codeml ancestral reconstruction (asr/rst)
  gbgc_genomic/      per-species NCOA1 genomic loci (~26 kb windows) + analysis scripts
  coactivator_scan/  nine-coactivator LXXLL scan inputs + output table
  foldx/             FoldX interaction-energy outputs for the NR-box–receptor complex
figures/             figure-generating scripts (build_gbgc_*.py)
```

## Dependencies

- Python 3 with `biopython`, `numpy`, `scipy`, `matplotlib`
- `mafft` (for the coactivator scan and the alignment-robustness check)
- HyPhy and PAML/codeml are required only to regenerate the selection results / `rst`
  from scratch; the committed outputs let you reproduce all figures and values without them.

## Reproduce

Run from the repository root. Scripts locate their inputs relative to this directory.

```bash
# All values reported in text (GC / percentile / CpG o/e; NR1-region composition;
# monotreme-stem W->S counts + Fisher's exact; NR1 site-level FUBAR/Contrast-FEL;
# platypus-echidna identities) — printed to stdout:
python analysis/gbgc_genomic/compute_reported_values.py

# Main figures:
python figures/build_gbgc_fig1.py          # motif + ancestral reconstruction
python figures/build_gbgc_fig2.py          # GC by clade + genomic tract
python figures/build_gbgc_fig3.py          # W->S bias + platypus-echidna retention
python figures/build_gbgc_fig4.py          # coactivator LXXLL scan

# Supplementary figures:
python figures/build_gbgc_figS_cpgisland.py   # CpG/GC enrichment across amniotes (S1)
python figures/build_gbgc_figS_codonpos.py    # GC bias by codon position (S2)
python figures/build_gbgc_figS_foldx.py       # FoldX NR-box interaction energy (S3)

# Verification (not a figure): NR1 call is robust to alignment method
python figures/build_gbgc_figS_alnrobust.py
```

## Re-fetching raw genomic loci (optional)

The genomic-locus FASTAs in `analysis/gbgc_genomic/` are ~26 kb windows fetched from public
NCBI/Ensembl accessions. To re-fetch them, set a contact email and run the fetch + QC scripts:

```bash
export ENTREZ_EMAIL="you@example.com"
python analysis/gbgc_genomic/fetch_reference_genomic.py
python analysis/gbgc_genomic/validate_genomic_fastas.py
```

## Source sequences

NCOA1 coding and genomic sequences are the public NCBI/Ensembl accessions for the 48 sampled
sarcopterygian species.
