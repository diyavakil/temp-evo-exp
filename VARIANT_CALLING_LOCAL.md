# Calling variants locally (fallback to the Plasmidsaurus website)

`scripts/run_variant_calling.sh` re-creates the steps Plasmidsaurus lists for its small-variant pipeline. Use it when you
need a reference the website cannot take (for example a newly sequenced ancestor), or if the website is unavailable.

> **It is an untested reconstruction.** Plasmidsaurus does not publish its exact Filtlong options, Clair3 model, Clair3
> options or filter thresholds, so those are placeholders here. Validate it (below) before using its output for
> anything you will publish.

## How the steps map

| Plasmidsaurus step | Local command | What is a guess |
|---|---|---|
| Filtlong v0.2.1, up to 100x | `filtlong --target_bases <genome size x 100>` | other Filtlong options are not stated |
| minimap2 v2.30 | `minimap2 -ax map-ont`, `samtools sort`, `samtools index` | preset `map-ont` is the standard choice |
| Clair3 v2.0.2, bacteria-tuned model | `run_clair3.sh --platform=ont --model_path=... --include_all_ctgs` | which model, and any extra options (set `CLAIR3_EXTRA`) |
| filter low-confidence and low-depth variants | `bcftools view -f PASS -i 'QUAL>=10 && FMT/DP>=20'` | the thresholds |
| SnpEff v5.4.0 | `snpEff build -genbank` then `snpEff` on the reference GenBank file | database built locally from your GBFF |
| structural variants | Sniffles2 (default) or assembly-vs-reference with minimap2 and paftools.js (`-V asm`) | the supplied tables are labelled sniffles2, but Plasmidsaurus's description names the assembly-based method |
| tables | `scripts/vcf_to_table.py` | no 'paralogous region' flag is produced |

## What you need

1. **Raw reads** (gzipped FASTQ). They are not in the standard results zip; request them from Plasmidsaurus.
2. **The reference FASTA** and its **GenBank file** (GBFF) with annotation. For ES lines, the NCBI ES114 assembly; for
   other lineages, the ancestor's assembly and its Bakta GBFF. Contig names in the FASTA and the GenBank LOCUS lines must
   match.
3. **A Clair3 bacteria-tuned model** that matches your basecalling. Ask Plasmidsaurus which model they used, and check
   Clair3's and Oxford Nanopore's model lists for the corresponding one. Their basecalling was Dorado v4.3
   super-accurate; I have not verified which Clair3 model corresponds to it.
4. The tools in `environment.yml`.

## Run

```bash
bash scripts/run_variant_calling.sh -s 18C_ESL1 \
  -r data/references/ES114_NCBI.fasta -g data/references/ES114_NCBI.gbff \
  -q data/reads/18C_ESL1.fastq.gz -m /path/to/clair3_model -o results/variant_calls/18C_ESL1 -t 8
```

Outputs in the output folder: the filtered reads, BAM, Clair3 folder, filtered and annotated VCFs, `tables/` with
`<sample> chr1.csv` etc., and `<sample>.versions.txt` recording tool versions and settings.

For a lineage with other contig names (for example the *V. jasicida* ET ancestor), pass `-M contigA=chr1,contigB=chr2,
contigC=pES100` so the table file names carry a replicon name that `subtract_variants.py` recognises.

## Validate before relying on it

Run the script on a sample for which you already have Plasmidsaurus tables (for example 18C ESL1 against the NCBI ES114
reference) and compare the two call sets with the same matching rules used for subtraction:

```bash
python scripts/subtract_variants.py \
  --a "results/variant_calls/18C_ESL1/tables/18C_ESL1 *.csv" --a-name local \
  --b "data/variant_tables/18C_ESL1 *.csv" --b-name plasmidsaurus --out results/validation
```

Most calls should be shared. Check the counts unique to each side, and adjust the Clair3 model, `-Q` and `-D` until
the local calls approximate Plasmidsaurus's. Note the differences in the write-up (for example the absent
paralogous-region flag and the structural-variant method).

## When to prefer the website

If the reference is already supported by the website, its tables are the simplest to cite and to reproduce. The local
script is mainly for new references, for transparency, and as a backup.
