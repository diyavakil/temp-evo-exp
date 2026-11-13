# Interpreting the k-mer results

## What a cell means

`containment_matrix.csv` has the header "query (row) in reference (column)". The value in row X, column Y is the share
of X's 31-mers that are also found in Y (1.0 = all, 0 = none).

The matrix is **directional**. If sample A contains two genomes and one of them is genome B, then B in A is close to 1.0
(everything in B is in A), but A in B is much lower (A has extra DNA). Large deletions and extra plasmids also make the
two directions differ. Read a pair as follows:

- both values high: the same genome;
- one high, one lower: one genome sits inside the other, which has extra sequence (a mixture, a plasmid or a deletion);
- both low: different genomes.

The values are a screen for **strain identity and large structural changes**. They cannot detect single-base
differences, so 0.9997 and 0.9998 cannot be told apart.

## Colour bands

| Containment | Colour | Meaning |
|---|---|---|
| 1.0 on the diagonal | gray | a genome compared with itself |
| 0.98 or higher | green | same strain |
| 0.90 to 0.98 | pale green | same strain with a large structural difference |
| 0.65 to 0.90 | yellow | unexpected middle range; review |
| 0.20 to 0.65 | orange | different strain of the same species, or a mixed assembly |
| below 0.20 | red | different species (the floor of about 0.002 comes from conserved sequence such as rRNA genes) |

Use four-decimal limits in Excel's "between" rules so ranges do not overlap, for example 0.98 to 1, 0.90 to 0.9799,
0.65 to 0.8999, 0.20 to 0.6499, 0 to 0.1999.

## Symmetric view and verdict formulas (Excel)

For a pair, judge by the **lower** of the two directions. On a new sheet, with the data in `Sheet1!B2:Y25`, fill this
across and down (adjust the range to your matrix size):

```
=MIN(Sheet1!B2, INDEX(Sheet1!$B$2:$Y$25, COLUMN()-1, ROW()-1))
```

A verdict for a single value in B2:

```
=IF(B2>=0.98,"same strain",IF(B2>=0.9,"same strain - review",IF(B2>=0.65,"review",IF(B2>=0.2,"different strain","different species"))))
```

## Figure caption (template)

> **Figure X. Genome identity between ancestral and evolved lines, measured by shared k-mers.** Each cell shows the
> fraction of the query genome's 31-mers (row) that are also found in the reference genome (column), where 1.0 means
> every k-mer is shared. Genomes were assembled from Oxford Nanopore reads, and k-mers were counted on both strands.
> Values are directional: a genome contained within a larger or mixed assembly scores high in one direction and low
> in the other. Cells are coloured by value: green (0.98 or higher) same strain; pale green (about 0.90 to 0.98) same
> strain with a large structural change; yellow (about 0.65 to 0.90) unexpected middle range (no cells fell here);
> orange (about 0.20 to 0.65) different strain of the same species or a mixed assembly; red (below 0.20) different
> species. Black outlines mark the lines expected to descend from each ancestor, based on line names.

## Cautions

- Low-coverage or fragmented assemblies can score lower than they should; check coverage and contig counts before
  reading a lower value as a biological difference.
- Species labels from Sourmash/GTDB do not separate strains of the same species (about 98% identity) and can give
  one label to an assembly that contains two genomes. The k-mer matrix and `kmer_identity.py contigs` can.
- Back up species-level statements in a paper with an alignment-based ANI (skani or FastANI) on the pairs that matter.
