#!/usr/bin/env python3
"""
convergence_estimate.py - rough chance estimate for 'the same gene was hit in more than one line'.

Model (deliberately simple): every gene is equally likely to receive a protein-changing mutation, each
line receives k_i such mutations, and gene length is ignored. The script reports how many genes you would
expect by chance to be hit in at least MIN_LINES lines, and the probability of seeing at least the observed
number of such genes (Poisson tail).

USAGE
  python convergence_estimate.py --genes 3800 --hits 7 6 7 --observed 3
  python convergence_estimate.py --genes 3800 --hits 7 6 7 --observed 3 --min-lines 3

This is a back-of-the-envelope screen, not a substitute for a formal test.
"""
import argparse
import math


def prob_hit_in_at_least(p_list, m):
    """P(a given gene is hit in at least m lines) given per-line hit probabilities (Poisson-binomial)."""
    dist = [1.0]
    for p in p_list:
        new = [0.0] * (len(dist) + 1)
        for j, d in enumerate(dist):
            new[j] += d * (1 - p)
            new[j + 1] += d * p
        dist = new
    return sum(dist[m:])


def poisson_tail(lam, k):
    """P(X >= k) for X ~ Poisson(lam)."""
    if k <= 0:
        return 1.0
    return 1.0 - sum(math.exp(-lam) * lam ** i / math.factorial(i) for i in range(k))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--genes", type=int, required=True, help="number of genes (CDS) in the genome")
    ap.add_argument("--hits", type=int, nargs="+", required=True, help="protein-changing hits in each line")
    ap.add_argument("--observed", type=int, required=True, help="genes observed hit in at least --min-lines lines")
    ap.add_argument("--min-lines", type=int, default=2)
    args = ap.parse_args()
    p = [k / args.genes for k in args.hits]
    q = prob_hit_in_at_least(p, args.min_lines)
    lam = args.genes * q
    print(f"lines: {len(args.hits)}; hits per line: {args.hits}; genes: {args.genes}")
    print(f"expected number of genes hit in >= {args.min_lines} lines by chance: {lam:.4f}")
    print(f"P(at least {args.observed} such genes by chance) = {poisson_tail(lam, args.observed):.2e}")


if __name__ == "__main__":
    main()
