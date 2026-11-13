# Methods and parameters

Parameters and decision rules used in the analysis, written so that each can be reproduced with the scripts in
this repository. Items in brackets need your own details.

## Sequencing and assembly (Plasmidsaurus)

Samples were sequenced and assembled by Plasmidsaurus: Dorado v4.3 super-accurate basecalling (default Q10 filter);
Filtlong v0.2.1 (reads under 300 bp removed, best 95% kept); genome size estimated with Lrge v0.2.1; assemblies built
with Autocycler from Flye v2.9.6+, Hifiasm and Plassembler v1.8.0+ and rotated with dnaapler; annotation with
Bakta v1.11; completeness with CheckM v1.2.2; species and plasmid identification with Sourmash against GTDB rs226,
RefSeq plasmids and phage sequences. [CONFIRM the Sourmash version and whether a polishing step was applied.]

## Identity verification (`scripts/kmer_identity.py`)

- **k-mers:** 31-mers, counted on both strands (a k-mer and its reverse complement count as one). Contigs of
  closed replicons are wrapped (ends joined); k-mers with ambiguous bases are discarded. Use `--linear` for
  draft assemblies of many linear contigs.
- **Subsampling:** one in eight k-mers is kept by hash value (FracMinHash-style) so that all-against-all comparison
  is fast. `--scale 1` gives exact values.
- **Containment:** for each ordered pair (query, reference), the fraction of the query's k-mers present in the
  reference. It is directional.
- **Approximate ANI:** containment raised to the power 1/31, times 100. Used only for orientation; decisions use
  containment.
- **Clustering:** single-linkage on the **lower** of the two directional values for each pair, threshold 0.90.
  A mixed assembly scores low in one direction against each parent, so it cannot join two genomes into one cluster.
- **Interpretation bands** (empirical, set from the observed distribution of values):

  | Containment | Meaning |
  |---|---|
  | 0.98 or higher | same strain |
  | 0.90 to 0.98 | same strain with a large structural difference (for example a deletion or an extra plasmid) |
  | 0.65 to 0.90 | unexpected; review |
  | 0.20 to 0.65 | different strain of the same species, or a mixed assembly |
  | below 0.20 | different species |

- **Decision rule:** a line proceeds to variant analysis only if it matches its expected ancestor in both directions
  (0.98 or higher). Lines matching another genome, or containing more than one, are flagged and excluded until the
  source is resolved.
- **Mixed samples:** `kmer_identity.py contigs` assigns each contig to the reference it matches.
- **Orthogonal confirmation** [PLANNED]: alignment-based ANI (skani or FastANI), whole-genome alignment (MUMmer
  `dnadiff`), a marker-gene check (gyrB, gapA or 16S) on the stocks, and mapping raw reads to candidate genomes.

## Variant calling

Calls were made by Plasmidsaurus against the uploaded reference: raw reads filtered with Filtlong to at most 100x,
aligned with minimap2 v2.30, small variants called with Clair3 v2.0.2 (bacteria-tuned model), low-confidence and
low-depth calls filtered, effects predicted with SnpEff v5.4.0. The pipeline is tuned for isolates, not for allele
fractions in mixed populations. [CONFIRM which caller produced the structural-variant tables (labelled sniffles2).]
`scripts/run_variant_calling.sh` is a local reconstruction of these steps and must be validated before use.

## Subtraction and triage (`scripts/subtract_variants.py`)

Two independently evolved lines at the same temperature are compared, and only calls found in one line are kept.

- **Matching:** SNPs must match on replicon, position and alternate base. Indels must be the same type on the same
  replicon within 10 bp (`--tol`). Structural variants must be the same type with start within 200 bp and end within
  10% of the size (at least 200 bp). Whole-replicon duplications caused by the circular reference are dropped.
- **Triage of what remains:**
  - *strong candidate*: VAF 0.9 or higher, depth 25 or higher, changes the protein (missense, nonsense, frameshift,
    in-frame change, start or stop lost; silent changes excluded), not within 200 bp of another call, and not in a
    region flagged as paralogous;
  - *likely artifact*: VAF below 0.7 and near a cluster of calls or in a flagged repeat; low-frequency structural
    variants;
  - *check*: everything else.
- **Dense clusters:** three or more calls chained within 200 bp are treated as one event.

## Mutator screen

For each line, count the unique SNPs and tabulate the substitution spectrum. A line is flagged as a probable mutator
when it carries many times more unique SNPs than its peers with one dominant class of change [SET: proposed rule, at
least 10 times the median and at least 90% in one class]. Confirm with a mutation-rate assay (for example a
rifampicin-resistance fluctuation test). Mutator lines are kept out of gene-level convergence counts.

## Convergence across lines

Genes with unique protein-changing calls in two or more independent non-mutator lines are tabulated and grouped by
pathway. `scripts/convergence_estimate.py` gives a rough chance estimate (equal mutation probability per gene,
gene length ignored). [PLANNED: a formal test.]

## Validation [PLANNED]

IGV inspection of candidate sites; PCR and Sanger sequencing of priority candidates and deletion junctions;
mutation-rate assay for suspected mutators; phenotypes against the ancestor (growth at each temperature,
luminescence, motility, biofilm, colonization).
