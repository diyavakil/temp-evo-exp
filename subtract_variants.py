#!/usr/bin/env python3
"""
subtract_variants.py - find variants that are UNIQUE to each of two evolved lines.

WHY
  If you call variants for every evolved line against the same reference (e.g. NCBI ES114),
  most calls are not new mutations. They are differences between your lab's ancestor and the
  reference sequence, plus recurring artifacts. Those show up in every line.
  Calls found in BOTH lines are therefore (almost certainly) not new. Calls found in only ONE
  line are the candidates for real evolved mutations.

WHAT IT NEEDS
  The variant tables (csv) you already have, one set per line, e.g.
     "18C_ESL1 chr1.csv", "18C_ESL1 chr2.csv", "18C_ESL1 pES100.csv"
  with the columns: source, POS, END, Variant Type, Variant Size, REF, ALT, DP, VAF,
  Affected Genes, Annotation, Impact, HGVS.p, Remark.
  The replicon (chr1 / chr2 / pES100) is taken from the file name.

USAGE  (put the patterns in quotes)
  python subtract_variants.py \
      --a "18C_ESL1 *.csv"  --a-name ESL1 \
      --b "VC 18C_ESL5 *.csv" --b-name ESL5 \
      --out results/

  Writes:
     <out>/unique_to_<A>.csv, <out>/unique_to_<B>.csv   (all unique calls, with a triage column)
     <out>/gene_overlap.csv                              (genes hit in BOTH lines = parallel evolution)

HOW CALLS ARE MATCHED
  SNP        same position and same ALT base.
  indel      same type (insertion / deletion) within --tol bp (default 10).
  structural same type, start within 200 bp, end within 10% of the size.

TRIAGE COLUMN
  strong candidate   present in nearly every read (VAF >= 0.9), depth >= 25, changes a protein
                     (not just 'upstream', and not a silent/synonymous change), and not sitting
                     next to a cluster of other calls
  check              unique, but missing one of those conditions (low VAF, noncoding, near a
                     cluster, flagged repeat region ...). Look at the reads before believing it.
  likely artifact    low VAF (<0.7) in a region where other calls cluster, or in a repeat
                     flagged by the pipeline

  These are rules of thumb to sort a list, NOT a verdict. Look at the reads (IGV) and confirm
  important ones by Sanger sequencing.
"""
import argparse
import glob
import os

import pandas as pd

CODING = ("missense_variant", "stop_gained", "frameshift_variant", "synonymous_variant",
          "disruptive_inframe", "conservative_inframe", "stop_lost", "start_lost",
          "inframe")


def load(pattern, label):
    parts = []
    files = sorted(glob.glob(pattern))
    if not files:
        raise SystemExit(f"no files match: {pattern}")
    for f in files:
        d = pd.read_csv(f)
        name = os.path.basename(f)
        d["rep"] = "pES100" if "pES100" in name else ("chr2" if "chr2" in name else "chr1")
        parts.append(d)
    df = pd.concat(parts, ignore_index=True)
    df["line"] = label
    return df


def is_match(q, pool, tol):
    p = pool[(pool.rep == q.rep) & (pool.source == q.source)]
    if q["source"] == "clair3":
        if q["Variant Type"] == "snp":
            m = p[(p.POS == q.POS) & (p.ALT == q.ALT)]
        else:
            m = p[(p["Variant Type"] == q["Variant Type"]) & ((p.POS - q.POS).abs() <= tol)]
    else:  # sniffles2 structural variants
        m = p[(p["Variant Type"] == q["Variant Type"]) & ((p.POS - q.POS).abs() <= 200)
              & ((p.END - q.END).abs() <= max(200, 0.1 * abs(q["Variant Size"])))]
    return len(m) > 0


def triage(row, df):
    ann = str(row.get("Annotation", ""))
    remark = str(row.get("Remark", ""))
    vaf, dp = row["VAF"], row["DP"]
    near = False
    if row["source"] == "clair3":
        p = df[(df.rep == row.rep) & (df.source == "clair3") & (df.index != row.name)]
        near = bool(((p.POS - row.POS).abs() <= 200).any())
    flagged_repeat = "paralogous" in remark
    coding = any(c in ann for c in CODING)
    if row["source"] == "sniffles2":
        return "check (structural variant: confirm in reads)" if vaf >= 0.5 else "likely artifact (low-frequency SV)"
    if vaf < 0.7 and (near or flagged_repeat):
        return "likely artifact"
    if ann.strip() == "synonymous_variant":
        return "check (silent change: protein unchanged)"
    if vaf >= 0.9 and dp >= 25 and coding and not near and not flagged_repeat:
        return "strong candidate"
    return "check"


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--a", required=True, help="glob for line A csv files (quote it)")
    ap.add_argument("--b", required=True, help="glob for line B csv files (quote it)")
    ap.add_argument("--a-name", default="A")
    ap.add_argument("--b-name", default="B")
    ap.add_argument("--tol", type=int, default=10, help="indel position tolerance in bp")
    ap.add_argument("--out", default="results")
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    A, B = load(args.a, args.a_name), load(args.b, args.b_name)
    for df, pool in ((A, B), (B, A)):
        df["shared"] = df.apply(lambda r: is_match(r, pool, args.tol), axis=1)

    keep = ["line", "rep", "source", "POS", "END", "Variant Type", "Variant Size", "REF", "ALT", "VAF", "DP",
            "Affected Genes", "Affected Products", "Annotation", "Impact", "HGVS.c", "HGVS.p", "Remark"]
    uniq = {}
    for df, name in ((A, args.a_name), (B, args.b_name)):
        c, s = df[df.source == "clair3"], df[df.source == "sniffles2"]
        print(f"{name}: {len(c)} small-variant calls -> {int(c.shared.sum())} shared with the other line, "
              f"{int((~c.shared).sum())} unique;  {len(s)} SVs -> {int(s.shared.sum())} shared, {int((~s.shared).sum())} unique")
        u = df[~df.shared].copy()
        # skip whole-replicon duplications caused by the circular reference
        u = u[~((u["Variant Type"] == "DUP") & (u.POS <= 1) & (u["Variant Size"] > 40000))]
        u["triage"] = u.apply(lambda r: triage(r, df), axis=1)
        u = u[keep + ["triage"]].sort_values(["rep", "POS"])
        u.to_csv(os.path.join(args.out, f"unique_to_{name}.csv"), index=False)
        uniq[name] = u
        print("   triage:", u.triage.value_counts().to_dict())

    # genes hit in both lines (parallel evolution candidates), using unique, non-artifact calls
    def genes(u):
        u = u[~u.triage.str.startswith("likely artifact")]
        g = {}
        for _, r in u.iterrows():
            if (isinstance(r["Affected Genes"], str) and any(c in str(r["Annotation"]) for c in CODING)
                    and str(r["Annotation"]).strip() != "synonymous_variant"):
                for gene in r["Affected Genes"].split("&"):
                    g.setdefault(gene, []).append(f"{r['Annotation']} {r['HGVS.p'] if isinstance(r['HGVS.p'], str) else ''}".strip())
        return g
    ga, gb = genes(uniq[args.a_name]), genes(uniq[args.b_name])
    both = sorted(set(ga) & set(gb))
    rows = [{"gene": g, f"{args.a_name}_change": "; ".join(ga[g]), f"{args.b_name}_change": "; ".join(gb[g])} for g in both]
    pd.DataFrame(rows).to_csv(os.path.join(args.out, "gene_overlap.csv"), index=False)
    print(f"\ngenes with a unique, protein-changing call in BOTH lines: {len(both)}")
    for r in rows:
        print("  ", r["gene"], "|", args.a_name + ":", r[f"{args.a_name}_change"], "|", args.b_name + ":", r[f"{args.b_name}_change"])


if __name__ == "__main__":
    main()
