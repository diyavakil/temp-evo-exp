# Genomic analysis of experimentally evolved *Vibrio* lines

Code for (1) checking that every evolved line still carries the genome it was supposed to, and (2) finding
mutations that arose during a serial-passage evolution experiment (lines evolved to generation 2000 at 18, 25
and 32 C, sequenced with Oxford Nanopore). [FILL IN: project title, authors, link to thesis or paper.]

> **Status: work in progress.** The identity-check and variant-filtering scripts were used to produce the
> results in the accompanying write-up. `scripts/run_variant_calling.sh` is an **untested reconstruction** of
> the Plasmidsaurus variant-calling pipeline (see [docs/VARIANT_CALLING_LOCAL.md](docs/VARIANT_CALLING_LOCAL.md)).

## Workflow

```
assemblies (.fasta)            raw reads (.fastq.gz)
      |                                |
      |  [gbff_to_fasta.py, only if you have a GenBank file instead of a FASTA]
      v                                |
 1. kmer_identity.py matrix            |   <- DECISION GATE: does each line match its expected ancestor?
      |  (lines that fail are excluded)|
      v                                v
 2. variant calling against the matching reference
      - Plasmidsaurus website (primary), or
      - scripts/run_variant_calling.sh (fallback; produces tables via vcf_to_table.py)
      |
      v
 3. subtract_variants.py   remove calls shared by two independent lines (ancestor-vs-reference differences,
      |                    recurring artifacts) and triage what remains
      v
 4. mutator screen, cross-line comparison, convergence_estimate.py, then validation (IGV, PCR, Sanger)
```

## Repository layout

```
scripts/
  kmer_identity.py          all-vs-all k-mer containment matrix + clustering; per-contig and difference modes
  gbff_to_fasta.py          GenBank -> FASTA with integrity checks (refuses damaged sequence)
  subtract_variants.py      find calls unique to each of two lines and triage them
  vcf_to_table.py           SnpEff-annotated VCF -> per-replicon CSV tables (same layout as Plasmidsaurus tables)
  run_variant_calling.sh    local variant-calling pipeline (UNTESTED reconstruction)
  convergence_estimate.py   rough chance estimate for 'same gene hit in several lines'
docs/
  METHODS.md                parameters and decision rules used in the analysis
  INTERPRETING_KMER_RESULTS.md   how to read the matrix, colour bands, Excel formulas, figure caption
  VARIANT_CALLING_LOCAL.md  how the local pipeline maps to the Plasmidsaurus steps, and how to validate it
tests/test_scripts.py       tests on synthetic data
data/                       NOT committed (see data/README.md); sample sheet template inside
results/                    output folder (not committed)
```

## Install

Python 3.8 or newer with `numpy` and `pandas` is enough for everything except the local variant-calling script:

```bash
pip install -r requirements.txt
python -m unittest discover -s tests -v      # 8 tests, a second or two
```

The local variant-calling script also needs Filtlong, minimap2, samtools, bcftools, Clair3, SnpEff and
(optionally) Sniffles2. `environment.yml` lists them for conda; Clair3 is easiest to run from its own container.

## 1. Identity check (k-mer containment)

```bash
mkdir -p assemblies results
cp path/to/*_1_reference.fasta assemblies/          # one FASTA per sample; the file name becomes the label
python scripts/kmer_identity.py matrix assemblies/*.fasta -o results/
```

Writes `assembly_stats.csv`, `containment_matrix.csv`, `ani_estimate_matrix.csv` and `clusters.csv`. Each cell of
the matrix is the share of the **row** genome's 31-mers that also occur in the **column** genome. Read
[docs/INTERPRETING_KMER_RESULTS.md](docs/INTERPRETING_KMER_RESULTS.md) before using the numbers.

- Add `--linear` when some files are draft assemblies made of many contigs (for example NCBI WGS records).
- Which contigs of one assembly come from which genome (finds mixed samples):
  `python scripts/kmer_identity.py contigs query.fasta refA.fasta refB.fasta`
- Which stretches of a reference are missing from a query (candidate deletions):
  `python scripts/kmer_identity.py diff reference.fasta query.fasta`
- GenBank file instead of FASTA: `python scripts/gbff_to_fasta.py sample.gbff --gff3 sample.gff3 -o assemblies/sample.fasta`

## 2. Variant calling

Use the Plasmidsaurus variant tables where you can. If you need to call variants yourself (a different
reference, a new ancestor, or if the website is unavailable), see
[docs/VARIANT_CALLING_LOCAL.md](docs/VARIANT_CALLING_LOCAL.md) and validate the local pipeline first.

## 3. Subtract and triage

```bash
python scripts/subtract_variants.py \
    --a "tables/18C_ESL1 *.csv" --a-name 18C_ESL1 \
    --b "tables/18C_ESL5 *.csv" --b-name 18C_ESL5 \
    --out results/18C
```

Compare **independent lines at the same temperature** (two lines evolved from the same ancestor). Calls found in
both are treated as differences between the lab stock and the reference, or recurring artifacts. The tables must be
named `<sample> chr1.csv`, `<sample> chr2.csv`, `<sample> pES100.csv` (the replicon is read from the file name).

## Known limitations

- The k-mer analysis is a screen for strain identity and large structural change. It cannot detect single-base
  differences; thresholds were set from the distribution of values in this data set.
- Subtraction against a second line also removes a real mutation that occurred at the same site in both lines.
  Calling against the lab's own sequenced ancestor is better and replaces the subtraction.
- Triage labels are rules of thumb for sorting a list, not verdicts. Confirm candidates in the aligned reads and
  by PCR/Sanger sequencing.
- Gene-level convergence estimates assume equal mutation probability per gene and ignore gene length.
