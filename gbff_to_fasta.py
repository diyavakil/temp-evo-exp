#!/usr/bin/env python3
"""
gbff_to_fasta.py - convert a GenBank file (.gbff / .gb / .gbk) to FASTA, with integrity checks.

Needs only Python 3.8+ (no Biopython).

WHY THIS EXISTS
  The k-mer script reads FASTA only. Bakta writes a .gbff (and a .gff3 with the sequence appended
  after a ##FASTA line). Files can be damaged in transfer, and a damaged sequence must never be
  converted silently: stray characters would simply be dropped and the k-mer results would be wrong.
  This script checks every replicon and refuses to write a FASTA unless every sequence is intact.

WHAT IT CHECKS (per replicon)
  - every sequence line is well formed (only nucleotide letters)
  - the sequence length equals the length declared in the file (LOCUS line, or ##sequence-region in a gff3)

USAGE
  python gbff_to_fasta.py sample.gbff -o sample.fasta
  python gbff_to_fasta.py sample.gbff --gff3 sample.gff3 -o sample.fasta
  python gbff_to_fasta.py sample.gff3 -o sample.fasta          # sequence is read from the ##FASTA section

  With --gff3, a replicon whose sequence is damaged in the .gbff is taken from the .gff3 instead, if
  it is intact there. This recovers a file in which the two copies are damaged in different places.
  Exit status is 0 when every replicon was written, 2 when any replicon could not be recovered
  (nothing is written in that case).

NOTE
  If Bakta produced the files, the assembly FASTA (.fna) with the same prefix holds the identical
  sequence and needs no conversion. Prefer it when you have it.
"""
import argparse
import re
import sys

NUC = re.compile(r"^[ACGTNRYKMSWBDHVacgtnrykmswbdhv]+$")
GB_LINE = re.compile(r"^\s*\d+(?: [A-Za-z]{1,10})+\s*$")


def parse_gbff(path):
    """Return {name: {declared, seq, bad_lines}} (file is read as latin-1 so it can never fail to decode)."""
    recs, name, in_origin = {}, None, False
    with open(path, encoding="latin-1") as fh:
        for line in fh:
            line = line.rstrip("\r\n")
            if line.startswith("LOCUS"):
                parts = line.split()
                name = parts[1]
                m = re.search(r"(\d+)\s+bp", line)
                recs[name] = {"declared": int(m.group(1)) if m else None, "seq": [], "bad_lines": 0}
                in_origin = False
            elif line.startswith("ORIGIN"):
                in_origin = True
            elif line.startswith("//"):
                in_origin = False
            elif in_origin and name and line.strip():
                if GB_LINE.match(line):
                    recs[name]["seq"].append(re.sub(r"[\s\d]", "", line))
                else:
                    recs[name]["bad_lines"] += 1
    return recs


def parse_gff3(path):
    """Return {name: {declared, seq, bad_lines}} from the ##FASTA section of a gff3 file."""
    declared, recs, name, in_fasta = {}, {}, None, False
    with open(path, encoding="latin-1") as fh:
        for line in fh:
            line = line.rstrip("\r\n")
            if not in_fasta:
                if line.startswith("##sequence-region"):
                    f = line.split()
                    if len(f) >= 4:
                        declared[f[1]] = int(f[3])
                elif line.startswith("##FASTA"):
                    in_fasta = True
                continue
            if line.startswith(">"):
                name = line[1:].split()[0]
                recs[name] = {"declared": declared.get(name), "seq": [], "bad_lines": 0}
            elif name and line.strip():
                if NUC.match(line.strip()):
                    recs[name]["seq"].append(line.strip())
                else:
                    recs[name]["bad_lines"] += 1
    return recs


def problems(rec):
    seq = "".join(rec["seq"])
    out = []
    if rec["bad_lines"]:
        out.append(f"{rec['bad_lines']} malformed sequence line(s)")
    if rec["declared"] is None:
        out.append("no declared length to check against")
    elif len(seq) != rec["declared"]:
        out.append(f"length {len(seq):,} != declared {rec['declared']:,}")
    if seq and not NUC.match(seq):
        out.append("non-nucleotide characters")
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("input", help=".gbff/.gb/.gbk file, or a .gff3 file with a ##FASTA section")
    ap.add_argument("--gff3", help="optional gff3 used to recover replicons that are damaged in the .gbff")
    ap.add_argument("-o", "--out", required=True, help="output FASTA")
    ap.add_argument("--width", type=int, default=80, help="line width of the output FASTA")
    args = ap.parse_args()

    primary = parse_gff3(args.input) if args.input.lower().endswith((".gff", ".gff3")) else parse_gbff(args.input)
    backup = parse_gff3(args.gff3) if args.gff3 else {}
    if not primary:
        sys.exit(f"no sequence records found in {args.input}")

    chosen, failed = {}, []
    for name, rec in primary.items():
        issues = problems(rec)
        if not issues:
            chosen[name] = ("".join(rec["seq"]).upper(), args.input)
            continue
        alt = backup.get(name)
        if alt is not None and not problems(alt):
            chosen[name] = ("".join(alt["seq"]).upper(), args.gff3)
            print(f"{name}: damaged in {args.input} ({'; '.join(issues)}); recovered from {args.gff3}")
        else:
            failed.append((name, issues, None if alt is None else problems(alt)))

    if failed:
        print("\nCould not recover every replicon; no FASTA written.", file=sys.stderr)
        for name, issues, alt_issues in failed:
            print(f"  {name}: {'; '.join(issues)}", file=sys.stderr)
            if alt_issues:
                print(f"     backup file also damaged: {'; '.join(alt_issues)}", file=sys.stderr)
        print("Get an undamaged copy of the file (or the assembly .fna) and try again.", file=sys.stderr)
        sys.exit(2)

    with open(args.out, "w") as out:
        for name, (seq, src) in chosen.items():
            out.write(f">{name}\n")
            for i in range(0, len(seq), args.width):
                out.write(seq[i:i + args.width] + "\n")
    print(f"wrote {len(chosen)} sequence(s) to {args.out}")
    for name, (seq, src) in chosen.items():
        print(f"  {name}: {len(seq):,} bp (from {src})")


if __name__ == "__main__":
    main()
