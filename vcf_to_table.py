#!/usr/bin/env python3
"""
vcf_to_table.py - turn a SnpEff-annotated VCF (e.g. from Clair3) into the per-replicon CSV tables used
by subtract_variants.py (the same layout as the Plasmidsaurus variant tables).

Needs only Python 3.8+ (standard library).

USAGE
  python vcf_to_table.py annotated.vcf.gz --sample 18C_ESL1 --out-dir tables/ \
         --rep-map NC_006840.2=chr1,NC_006841.2=chr2,NC_006842.1=pES100

  Writes one file per replicon, named "<sample> <rep>.csv" (for example "18C_ESL1 chr1.csv"), because
  subtract_variants.py reads the replicon (chr1 / chr2 / pES100) from the file name.

COLUMNS WRITTEN
  source, POS, END, Variant Type (snp / mnp / insertion / deletion), Variant Size, REF, ALT, QUAL, FILTER,
  DP, VAF, GT, Affected Genes, Affected Products, Annotation, Impact, HGVS.c, HGVS.p, Remark

NOT REPRODUCED
  The Plasmidsaurus 'Remark' column (paralogous-region flag) and 'Affected Products' are left empty.
"""
import argparse
import csv
import gzip
import os
import sys

COLUMNS = ["source", "POS", "END", "Variant Type", "Variant Size", "REF", "ALT", "QUAL", "FILTER", "DP", "VAF",
           "GT", "Affected Genes", "Affected Products", "Annotation", "Impact", "HGVS.c", "HGVS.p", "Remark"]
DEFAULT_MAP = "NC_006840.2=chr1,NC_006841.2=chr2,NC_006842.1=pES100"


def open_text(path):
    return gzip.open(path, "rt", encoding="latin-1") if path.endswith(".gz") else open(path, encoding="latin-1")


def variant_type(ref, alt):
    if len(ref) == 1 and len(alt) == 1:
        return "snp"
    if len(ref) == len(alt):
        return "mnp"
    return "insertion" if len(alt) > len(ref) else "deletion"


def parse_ann(info):
    for item in info.split(";"):
        if item.startswith("ANN="):
            entries = [e.split("|") for e in item[4:].split(",")]
            first = entries[0] + [""] * (11 - len(entries[0]))
            genes = []
            for e in entries:
                if len(e) > 3 and e[3] and e[3] not in genes:
                    genes.append(e[3])
            return first[1], first[2], "&".join(genes), first[9], first[10]
    return "", "", "", "", ""


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("vcf", help="annotated VCF (.vcf or .vcf.gz)")
    ap.add_argument("--sample", required=True, help="name used in the output file names")
    ap.add_argument("--out-dir", default="tables")
    ap.add_argument("--rep-map", default=DEFAULT_MAP,
                    help=f"contig=replicon pairs, comma separated (default: {DEFAULT_MAP}). Replicon names must be "
                         "chr1, chr2 or pES100 for subtract_variants.py to recognise them")
    ap.add_argument("--source", default="clair3", help="value for the 'source' column (default clair3)")
    ap.add_argument("--min-qual", type=float, default=0.0)
    ap.add_argument("--min-dp", type=int, default=0)
    ap.add_argument("--pass-only", action="store_true", help="keep only FILTER=PASS (or '.') records")
    args = ap.parse_args()

    repmap = dict(p.split("=") for p in args.rep_map.split(",") if p)
    os.makedirs(args.out_dir, exist_ok=True)
    rows, skipped = {}, 0
    with open_text(args.vcf) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 8:
                continue
            chrom, pos, _id, ref, alt, qual, flt, info = f[:8]
            rep = repmap.get(chrom)
            if rep is None:
                skipped += 1
                continue
            if args.pass_only and flt not in ("PASS", "."):
                continue
            q = float(qual) if qual not in (".", "") else 0.0
            if q < args.min_qual:
                continue
            alt = alt.split(",")[0]
            gt = dp = vaf = ""
            if len(f) >= 10:
                keys = f[8].split(":")
                vals = dict(zip(keys, f[9].split(":")))
                gt = vals.get("GT", "")
                dp = vals.get("DP", "")
                if "AF" in vals and vals["AF"] not in (".", ""):
                    vaf = vals["AF"].split(",")[0]
                elif "AD" in vals and dp:
                    ad = vals["AD"].split(",")
                    if len(ad) >= 2 and int(dp) > 0:
                        vaf = f"{int(ad[1]) / int(dp):.4f}"
            if args.min_dp and (not dp or int(dp) < args.min_dp):
                continue
            annot, impact, genes, hgvs_c, hgvs_p = parse_ann(info)
            p = int(pos)
            rows.setdefault(rep, []).append({
                "source": args.source, "POS": p, "END": p + len(ref) - 1, "Variant Type": variant_type(ref, alt),
                "Variant Size": len(alt) - len(ref), "REF": ref, "ALT": alt, "QUAL": q, "FILTER": flt, "DP": dp,
                "VAF": vaf, "GT": gt, "Affected Genes": genes, "Affected Products": "", "Annotation": annot,
                "Impact": impact, "HGVS.c": hgvs_c, "HGVS.p": hgvs_p, "Remark": ""})
    for rep, recs in rows.items():
        path = os.path.join(args.out_dir, f"{args.sample} {rep}.csv")
        with open(path, "w", newline="") as out:
            w = csv.DictWriter(out, fieldnames=COLUMNS)
            w.writeheader()
            w.writerows(sorted(recs, key=lambda r: r["POS"]))
        print(f"{path}: {len(recs)} variants")
    if skipped:
        print(f"{skipped} record(s) on contigs not in --rep-map were skipped", file=sys.stderr)


if __name__ == "__main__":
    main()
