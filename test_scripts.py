"""Self-contained tests for the analysis scripts (synthetic data only).

Run from the repository root:   python -m unittest discover -s tests -v
"""
import csv
import glob
import os
import random
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
SCRIPTS = ROOT / "scripts"
sys.path.insert(0, str(SCRIPTS))


def run(*args, check=True):
    return subprocess.run([sys.executable, *map(str, args)], capture_output=True, text=True, check=check)


def random_genome(n, seed):
    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(n))


def revcomp(s):
    return s.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def write_fasta(path, records):
    with open(path, "w") as fh:
        for name, seq in records.items():
            fh.write(f">{name}\n")
            for i in range(0, len(seq), 80):
                fh.write(seq[i:i + 80] + "\n")


def read_matrix(path):
    with open(path) as fh:
        rows = list(csv.reader(fh))
    cols = rows[0][1:]
    return {r[0]: dict(zip(cols, map(float, r[1:]))) for r in rows[1:]}


class KmerIdentityTests(unittest.TestCase):
    def test_matrix_and_clusters(self):
        with tempfile.TemporaryDirectory() as d:
            a = random_genome(60000, 1)
            b = a[:20000] + a[23000:]          # A with a 3 kb deletion
            c = revcomp(a)                      # same genome, other strand
            e = random_genome(60000, 2)         # unrelated
            for name, seq in [("A", a), ("B", b), ("C", c), ("D", e)]:
                write_fasta(f"{d}/{name}.fasta", {"chr": seq})
            run(SCRIPTS / "kmer_identity.py", "matrix", *[f"{d}/{n}.fasta" for n in "ABCD"],
                "--scale", 1, "-o", f"{d}/out")
            m = read_matrix(f"{d}/out/containment_matrix.csv")
            self.assertAlmostEqual(m["A"]["C"], 1.0, places=4)   # strand independent
            self.assertGreater(m["B"]["A"], 0.99)                # deletion: B sits inside A ...
            self.assertAlmostEqual(m["A"]["B"], 0.95, delta=0.01)  # ... but A is not fully inside B
            self.assertLess(m["A"]["D"], 0.01)                   # unrelated genomes
            with open(f"{d}/out/clusters.csv") as fh:
                clusters = {}
                for r in list(csv.reader(fh))[1:]:
                    clusters.setdefault(r[0], []).append(r[1])
            groups = sorted(sorted(v) for v in clusters.values())
            self.assertEqual(groups, [["A", "B", "C"], ["D"]])

    def test_mixed_assembly_does_not_join_clusters(self):
        with tempfile.TemporaryDirectory() as d:
            a, e = random_genome(50000, 3), random_genome(50000, 4)
            write_fasta(f"{d}/A.fasta", {"chr": a})
            write_fasta(f"{d}/E.fasta", {"chr": e})
            write_fasta(f"{d}/MIX.fasta", {"c1": a, "c2": e})     # contains both genomes
            run(SCRIPTS / "kmer_identity.py", "matrix", f"{d}/A.fasta", f"{d}/E.fasta", f"{d}/MIX.fasta",
                "--scale", 1, "-o", f"{d}/out")
            m = read_matrix(f"{d}/out/containment_matrix.csv")
            self.assertGreater(m["A"]["MIX"], 0.99)
            self.assertAlmostEqual(m["MIX"]["A"], 0.5, delta=0.02)
            with open(f"{d}/out/clusters.csv") as fh:
                n_clusters = len({r[0] for r in list(csv.reader(fh))[1:]})
            self.assertEqual(n_clusters, 3)                       # the mixture links nothing

    def test_diff_finds_deletion(self):
        with tempfile.TemporaryDirectory() as d:
            a = random_genome(80000, 5)
            b = a[:30000] + a[40000:]                              # 10 kb deletion
            write_fasta(f"{d}/ref.fasta", {"chr": a})
            write_fasta(f"{d}/qry.fasta", {"chr": b})
            out = run(SCRIPTS / "kmer_identity.py", "diff", f"{d}/ref.fasta", f"{d}/qry.fasta").stdout
            self.assertIn("candidate deletion", out)
            self.assertRegex(out, r"29,9\d\d-40,0\d\d")           # junction k-mers start ~30 bp before 30 kb


def make_gbff(name, seq, corrupt=False):
    lines = [f"LOCUS       {name:<16}{len(seq):>9} bp    DNA     linear   UNK 01-JAN-2026", "ORIGIN"]
    for i in range(0, len(seq), 60):
        chunk = seq[i:i + 60].lower()
        groups = " ".join(chunk[j:j + 10] for j in range(0, len(chunk), 10))
        lines.append(f"{i + 1:>9} {groups}")
    if corrupt:
        lines[4] = "        181 \x93\xa0\x02acgt\x7f\xa0 ??? zz"
    lines.append("//")
    return "\n".join(lines) + "\n"


def make_gff3(seqs, corrupt=()):
    out = ["##gff-version 3"]
    out += [f"##sequence-region {n} 1 {len(s)}" for n, s in seqs.items()]
    out.append("##FASTA")
    for n, s in seqs.items():
        out.append(f">{n}")
        body = [s[i:i + 60] for i in range(0, len(s), 60)]
        if n in corrupt:
            body[2] = "AC\x93\xa0GT??"
        out += body
    return "\n".join(out) + "\n"


class GbffToFastaTests(unittest.TestCase):
    def setUp(self):
        self.s1, self.s2 = random_genome(600, 11), random_genome(450, 12)

    def test_clean_file(self):
        with tempfile.TemporaryDirectory() as d:
            Path(f"{d}/x.gbff").write_text(make_gbff("r1", self.s1) + make_gbff("r2", self.s2))
            run(SCRIPTS / "gbff_to_fasta.py", f"{d}/x.gbff", "-o", f"{d}/x.fasta")
            text = Path(f"{d}/x.fasta").read_text().split(">")[1:]
            got = {t.split("\n")[0]: "".join(t.split("\n")[1:]) for t in text}
            self.assertEqual(got, {"r1": self.s1, "r2": self.s2})

    def test_recovers_damaged_replicon_from_gff3(self):
        with tempfile.TemporaryDirectory() as d:
            Path(f"{d}/x.gbff").write_text(make_gbff("r1", self.s1) + make_gbff("r2", self.s2, corrupt=True))
            Path(f"{d}/x.gff3").write_text(make_gff3({"r1": self.s1, "r2": self.s2}, corrupt=("r1",)))
            run(SCRIPTS / "gbff_to_fasta.py", f"{d}/x.gbff", "--gff3", f"{d}/x.gff3", "-o", f"{d}/x.fasta")
            text = Path(f"{d}/x.fasta").read_text().split(">")[1:]
            got = {t.split("\n")[0]: "".join(t.split("\n")[1:]) for t in text}
            self.assertEqual(got, {"r1": self.s1, "r2": self.s2})

    def test_refuses_when_unrecoverable(self):
        with tempfile.TemporaryDirectory() as d:
            Path(f"{d}/x.gbff").write_text(make_gbff("r1", self.s1, corrupt=True))
            result = run(SCRIPTS / "gbff_to_fasta.py", f"{d}/x.gbff", "-o", f"{d}/x.fasta", check=False)
            self.assertEqual(result.returncode, 2)
            self.assertFalse(os.path.exists(f"{d}/x.fasta"))      # nothing written


VCF = """##fileformat=VCFv4.2
##contig=<ID=NC_006840.2>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS
{rows}
"""


def vcf_row(chrom, pos, ref, alt, gene="geneA", ann="missense_variant", impact="MODERATE", af="1.0", dp=50):
    info = f"ANN={alt}|{ann}|{impact}|{gene}|{gene}|transcript|T1|protein_coding|1/1|c.1A>T|p.Met1Leu||||||"
    return f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t18\tPASS\t{info}\tGT:DP:AF\t1/1:{dp}:{af}"


class VariantTableTests(unittest.TestCase):
    def test_vcf_to_table_and_subtraction(self):
        with tempfile.TemporaryDirectory() as d:
            shared = vcf_row("NC_006840.2", 1000, "A", "G", gene="shared")
            a_only = vcf_row("NC_006840.2", 2000, "C", "T", gene="onlyA")
            b_only = vcf_row("NC_006840.2", 3000, "GA", "G", gene="onlyB", ann="frameshift_variant", impact="HIGH")
            Path(f"{d}/a.vcf").write_text(VCF.format(rows="\n".join([shared, a_only])))
            Path(f"{d}/b.vcf").write_text(VCF.format(rows="\n".join([shared, b_only])))
            for s in "ab":
                run(SCRIPTS / "vcf_to_table.py", f"{d}/{s}.vcf", "--sample", f"line{s.upper()}", "--out-dir", f"{d}/t")
            with open(f"{d}/t/lineB chr1.csv") as fh:
                rows = list(csv.DictReader(fh))
            self.assertEqual([r["Variant Type"] for r in rows], ["snp", "deletion"])
            self.assertEqual(rows[1]["Variant Size"], "-1")
            self.assertEqual(rows[0]["Affected Genes"], "shared")
            self.assertEqual(rows[0]["VAF"], "1.0")
            run(SCRIPTS / "subtract_variants.py", "--a", f"{d}/t/lineA *.csv", "--a-name", "A",
                "--b", f"{d}/t/lineB *.csv", "--b-name", "B", "--out", f"{d}/res")
            with open(f"{d}/res/unique_to_A.csv") as fa, open(f"{d}/res/unique_to_B.csv") as fb:
                ua, ub = list(csv.DictReader(fa)), list(csv.DictReader(fb))
            self.assertEqual([r["Affected Genes"] for r in ua], ["onlyA"])
            self.assertEqual([r["Affected Genes"] for r in ub], ["onlyB"])


class ConvergenceTests(unittest.TestCase):
    def test_expected_value(self):
        import convergence_estimate as ce
        p = [7 / 3800, 6 / 3800, 7 / 3800]
        lam = 3800 * ce.prob_hit_in_at_least(p, 2)
        self.assertAlmostEqual(lam, 0.035, places=3)
        self.assertLess(ce.poisson_tail(lam, 3), 1e-5)


if __name__ == "__main__":
    unittest.main()
