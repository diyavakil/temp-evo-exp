#!/usr/bin/env python3
"""
kmer_identity.py - quick strain-identity and structural checks between genome assemblies.

Needs only Python 3.8+ and numpy  (pip install numpy).  No databases, no internet.

WHY THIS EXISTS
  Before calling SNPs between an evolved line and its ancestor you want to know that they
  really are the same strain.  Same-strain assemblies share >99% of their 31-mers.
  Different strains of the same species (~98% ANI) share only ~55-60%.  Different species
  share almost none.  This script measures that directly.

THREE COMMANDS
  1) matrix   all-vs-all comparison + clustering of many assemblies
        python kmer_identity.py matrix  assemblies/*.fasta  -o results/

        writes:  assembly_stats.csv         contigs, lengths, GC, circular flags
                 containment_matrix.csv     row = query, column = reference; value = fraction of the
                                            row genome's k-mers that are found in the column genome
                 ani_estimate_matrix.csv    same matrix converted to approximate ANI (%)
                 clusters.csv               which genomes group together as "same strain"

  2) contigs  which reference does each contig of a query come from?  (finds mixed samples)
        python kmer_identity.py contigs  query.fasta  refA.fasta refB.fasta ...

  3) diff     where is query missing sequence that the reference has?  (deletions / divergence)
        python kmer_identity.py diff  reference.fasta  query.fasta  [--window 100000] [--min-run 5000]
        Reports, per window along the reference, the % of reference k-mers absent from the query,
        and lists long runs of missing sequence (candidate deletions, in reference coordinates).

HOW TO READ THE NUMBERS
  containment(A in B) = share of A's k-mers present in B.
     ~0.99+   same strain (differences = a handful of SNPs / indels)
     ~0.55-0.6 different strains of one species (~98% ANI)
     ~0.0     unrelated
  approximate ANI = containment ** (1/k).   Rough, but good for sorting same / different strain.
  A large deletion in one genome makes containment asymmetric (the genome that LOST sequence is
  still fully contained in the other, but not vice versa), so clustering uses the MINIMUM of the
  two directions.  A mixed sample (two genomes in one assembly) shows ~0.4-0.6 in both of its
  parents: use the `contigs` command on it.

LIMITS
  This is a k-mer screen, not an alignment.  It will not tell you which base changed.  For SNPs,
  map reads to the ancestor assembly and call variants from the pileup.
"""
import argparse
import csv
import os
import re
import sys

import numpy as np

K = 31  # k-mer size (odd, <= 31 so a k-mer fits in 62 bits)

_TBL = np.full(256, 255, dtype=np.uint8)
for _i, _c in enumerate("ACGT"):
    _TBL[ord(_c)] = _i


# ----------------------------------------------------------------------------------------
# I/O helpers
# ----------------------------------------------------------------------------------------
def read_fasta(path):
    seqs, cur = {}, None
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                cur = line[1:]
                seqs[cur] = []
            elif cur is not None:
                seqs[cur].append(line)
    return {name: "".join(parts).upper() for name, parts in seqs.items()}


def short_label(path):
    b = os.path.basename(path)
    b = re.sub(r"\.(fasta|fa|fna)$", "", b)
    b = re.sub(r"(_1_reference(_1)?)$", "", b)
    b = re.sub(r"_2000gn$", "", b)
    return b


def contig_id(header):
    return header.split()[0]


def is_circular(header):
    return "circular=true" in header


# ----------------------------------------------------------------------------------------
# k-mer machinery (canonical k-mers, i.e. a k-mer and its reverse complement count as one)
# ----------------------------------------------------------------------------------------
def canonical_kmers(seq, circular=True):
    """Array of canonical k-mers (uint64) for one sequence.  Positions follow the sequence."""
    if circular:
        seq = seq + seq[: K - 1]
    a = _TBL[np.frombuffer(seq.encode(), dtype=np.uint8)].astype(np.uint64)
    n = len(a) - K + 1
    if n <= 0:
        return np.zeros(0, dtype=np.uint64)
    fw = np.zeros(n, dtype=np.uint64)
    rc = np.zeros(n, dtype=np.uint64)
    for i in range(K):
        fw = (fw << np.uint64(2)) | a[i:i + n]
        rc = rc | ((np.uint64(3) - a[i:i + n]) << np.uint64(2 * i))
    kmers = np.minimum(fw, rc)
    # k-mers containing N / non-ACGT bases were encoded from 255 -> discard garbage by masking
    bad = (a == np.uint64(255))
    if bad.any():
        cs = np.concatenate([[0], np.cumsum(bad)])
        has_bad = (cs[K:K + n] - cs[:n]) > 0
        kmers = kmers[~has_bad]
    return kmers


def subsample(kmers, scale):
    """FracMinHash-style subsample: keep ~1/scale of k-mers, deterministically by hash value."""
    if scale <= 1:
        return kmers
    h = kmers * np.uint64(0x9E3779B97F4A7C15)
    return kmers[((h >> np.uint64(8)) % np.uint64(scale)) == 0]


def contain(a_sorted_unique, b_sorted_unique):
    """Fraction of a's k-mers that appear in b (both sorted unique uint64 arrays)."""
    if len(a_sorted_unique) == 0 or len(b_sorted_unique) == 0:
        return float("nan")
    idx = np.searchsorted(b_sorted_unique, a_sorted_unique)
    idx[idx == len(b_sorted_unique)] = len(b_sorted_unique) - 1
    return float((b_sorted_unique[idx] == a_sorted_unique).mean())


def genome_kmers(path, scale, circular=True):
    d = read_fasta(path)
    allk = np.concatenate([canonical_kmers(s, circular=circular) for s in d.values()])
    return np.unique(subsample(allk, scale)), d


# ----------------------------------------------------------------------------------------
# command: matrix
# ----------------------------------------------------------------------------------------
def cmd_matrix(args):
    paths = []
    seen = set()
    for p in args.fasta:
        rp = os.path.realpath(p)
        if rp not in seen:
            seen.add(rp)
            paths.append(p)
    os.makedirs(args.out, exist_ok=True)

    names, sets, stats = [], [], []
    for p in paths:
        ks, d = genome_kmers(p, args.scale, circular=not args.linear)
        name = short_label(p)
        names.append(name)
        sets.append(ks)
        total = sum(len(s) for s in d.values())
        gc = sum(s.count("G") + s.count("C") for s in d.values()) / total * 100
        lens = sorted((len(s) for s in d.values()), reverse=True)
        flags = "".join("c" if is_circular(h) else "L" for h in d)
        stats.append([name, len(d), total, round(gc, 2), ";".join(map(str, lens)), flags])
        print(f"  loaded {name}: {len(d)} contigs, {total:,} bp", file=sys.stderr)

    n = len(names)
    C = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            C[i, j] = 1.0 if i == j else contain(sets[i], sets[j])
    ANI = np.where(C > 0, C ** (1.0 / K) * 100, 0.0)

    with open(os.path.join(args.out, "assembly_stats.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["assembly", "n_contigs", "total_bp", "GC_percent",
                    "contig_lengths_bp", "circular_flags (c=circular, L=linear, in contig order)"])
        w.writerows(stats)

    for fname, M, fmt in [("containment_matrix.csv", C, "{:.4f}"),
                          ("ani_estimate_matrix.csv", ANI, "{:.2f}")]:
        with open(os.path.join(args.out, fname), "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["query (row) in reference (column)"] + names)
            for i in range(n):
                w.writerow([names[i]] + [fmt.format(M[i, j]) for j in range(n)])

    # single-linkage clustering on the symmetric (minimum) containment
    S = np.minimum(C, C.T)
    parent = list(range(n))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for i in range(n):
        for j in range(i + 1, n):
            if S[i, j] >= args.threshold:
                parent[find(i)] = find(j)
    groups = {}
    for i in range(n):
        groups.setdefault(find(i), []).append(i)
    ordered = sorted(groups.values(), key=lambda g: (-len(g), names[g[0]]))

    with open(os.path.join(args.out, "clusters.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["cluster", "assembly", "lowest_symmetric_containment_with_cluster_mates"])
        for ci, g in enumerate(ordered, 1):
            for i in sorted(g, key=lambda i: names[i]):
                mates = [S[i, j] for j in g if j != i]
                w.writerow([ci, names[i], f"{min(mates):.3f}" if mates else "singleton"])

    print(f"\nWrote results to {args.out}/", file=sys.stderr)
    for ci, g in enumerate(ordered, 1):
        print(f"cluster {ci}: " + ", ".join(sorted(names[i] for i in g)))


# ----------------------------------------------------------------------------------------
# command: contigs
# ----------------------------------------------------------------------------------------
def cmd_contigs(args):
    q = read_fasta(args.query)
    refs = {}
    for r in args.refs:
        ks, _ = genome_kmers(r, 1)
        refs[short_label(r)] = ks
    print(f"{'contig':<14}{'length':>10}  " + "  ".join(f"{n[:18]:>18}" for n in refs))
    for h, s in q.items():
        k = np.unique(canonical_kmers(s))
        row = "  ".join(f"{contain(k, ks):>18.3f}" for ks in refs.values())
        print(f"{contig_id(h):<14}{len(s):>10,}  {row}")
    print("\n(values = fraction of the contig's k-mers found in each reference; "
          "~1.0 means the contig belongs to that genome)")


# ----------------------------------------------------------------------------------------
# command: diff
# ----------------------------------------------------------------------------------------
def cmd_diff(args):
    ref = read_fasta(args.reference)
    qry = read_fasta(args.query)
    qk = np.unique(np.concatenate([canonical_kmers(s) for s in qry.values()]))
    print(f"reference: {args.reference}\nquery:     {args.query}")
    print("k-mers of the REFERENCE that are absent from the QUERY (per window along each contig)\n")
    for h, s in ref.items():
        km = canonical_kmers(s)
        idx = np.searchsorted(qk, km)
        idx[idx == len(qk)] = len(qk) - 1
        missing = qk[idx] != km
        print(f"{contig_id(h)}  length={len(s):,}  missing={missing.mean() * 100:.2f}%")
        wins = []
        for start in range(0, len(km), args.window):
            wins.append(f"{start // 1000}k:{missing[start:start + args.window].mean() * 100:.1f}%")
        print("   " + "  ".join(wins))
        pos = np.where(missing)[0]
        if len(pos):
            runs, st, prev = [], pos[0], pos[0]
            for p in pos[1:]:
                if p - prev > args.gap:
                    runs.append((st, prev))
                    st = p
                prev = p
            runs.append((st, prev))
            for a, b in runs:
                if b - a >= args.min_run:
                    print(f"   candidate deletion in query: {a:,}-{b + K - 1:,} "
                          f"(~{(b - a + K) / 1000:.1f} kb, reference coordinates)")
        print()


# ----------------------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)

    m = sub.add_parser("matrix", help="all-vs-all identity matrix and clusters")
    m.add_argument("fasta", nargs="+")
    m.add_argument("-o", "--out", default="results")
    m.add_argument("--scale", type=int, default=8,
                   help="keep 1/scale of k-mers (default 8: fast, still accurate for this purpose; 1 = exact)")
    m.add_argument("--linear", action="store_true",
                   help="treat contigs as linear (do not wrap the ends). Use for draft assemblies made of many "
                        "contigs (e.g. NCBI WGS records); default is to wrap, which suits closed circular replicons")
    m.add_argument("--threshold", type=float, default=0.90,
                   help="symmetric containment needed to call two genomes the same strain (default 0.90)")
    m.set_defaults(func=cmd_matrix)

    c = sub.add_parser("contigs", help="assign each contig of a query to the reference it matches")
    c.add_argument("query")
    c.add_argument("refs", nargs="+")
    c.set_defaults(func=cmd_contigs)

    d = sub.add_parser("diff", help="find sequence present in reference but missing from query")
    d.add_argument("reference")
    d.add_argument("query")
    d.add_argument("--window", type=int, default=100000)
    d.add_argument("--min-run", type=int, default=5000, help="smallest deletion to report (bp)")
    d.add_argument("--gap", type=int, default=2000, help="merge missing stretches closer than this (bp)")
    d.set_defaults(func=cmd_diff)

    args = ap.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
