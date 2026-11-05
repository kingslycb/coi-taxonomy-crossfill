#!/usr/bin/env python3
"""
Convert the raw output of each classifier into one table per library with the same columns:
ASV, kingdom ... species, and the method's support values. Unassigned ranks are 'NA'; once a rank fails its
support threshold, all finer ranks are unassigned too.

  blast   BLASTn tabular output (-outfmt "6 qseqid sseqid pident length qlen slen evalue bitscore"), top hit per ASV;
          the hit's lineage (ref_lineages.tsv) is kept down to the finest rank whose identity threshold is met.
  sintax  vsearch --sintax tab-separated output; ranks with bootstrap support below the cut-off are unassigned.
  rdp     RDP Classifier 'allrank' output; ranks with confidence below the cut-off are unassigned.

Usage
  classification_tables.py blast  --asvs A.fa --raw A.blast.tsv --lineages ref_lineages.tsv --identity 98,95,90,85,80,75 --out A_BLAST.tsv
  classification_tables.py sintax --asvs A.fa --raw A.sintax.tsv --cutoff 0.75 --out A_SINTAX.tsv
  classification_tables.py rdp    --asvs A.fa --raw A.rdp.tsv    --cutoff 0.75 --out A_RDP.tsv
"""
import argparse
import re

from taxutils import RANKS, SINTAX_CODE, read_fasta

NA = "NA"


def asv_ids(path):
    return [header.split()[0] for header, _ in read_fasta(path)]


def truncate(names, keep):
    """names, keep: lists over RANKS. A rank is assigned only if its support and the support of every coarser rank
    pass; a rank missing from the reference lineage stays NA without blocking the finer ranks."""
    out, ok = [], True
    for name, k in zip(names, keep):
        ok = ok and k
        out.append(name if ok and name else NA)
    return out


def write_table(path, rows, extra):
    with open(path, "w") as f:
        f.write("\t".join(["ASV"] + RANKS + extra) + "\n")
        for row in rows:
            f.write("\t".join(str(x) for x in row) + "\n")


def blast(args):
    lineage = {}
    with open(args.lineages) as f:
        cols = f.readline().rstrip("\n").split("\t")
        for line in f:
            p = dict(zip(cols, line.rstrip("\n").split("\t")))
            if p["database"] == "MIDORI2":
                lineage[p["ref_id"]] = p
    thresholds = dict(zip(["species", "genus", "family", "order", "class", "phylum"],
                          (float(x) for x in args.identity.split(","))))
    best = {}
    with open(args.raw) as f:
        for line in f:
            q, s, pident, length, qlen, slen, evalue, bitscore = line.rstrip("\n").split("\t")[:8]
            key = (float(bitscore), float(pident))
            if q not in best or key > best[q][0]:
                best[q] = (key, s, float(pident), evalue)
    rows = []
    for asv in asv_ids(args.asvs):
        if asv not in best:
            rows.append([asv] + [NA] * len(RANKS) + [NA, NA, NA, NA])
            continue
        _, ref, pident, evalue = best[asv]
        lin = lineage[ref]
        keep = [pident >= thresholds.get(r, 0.0) for r in RANKS]
        rows.append([asv] + truncate([lin[r] for r in RANKS], keep) + [f"{pident:.2f}", evalue, ref,
                                                                       lin["source_database"]])
    write_table(args.out, rows, ["percent_identity", "evalue", "reference_id", "reference_source"])


def sintax(args):
    code_rank = {v: k for k, v in SINTAX_CODE.items()}
    parsed = {}
    with open(args.raw) as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            names, conf = {}, {}
            for part in (p[1].split(",") if len(p) > 1 and p[1] else []):
                m = re.match(r"([a-z]):(.*)\(([\d.]+)\)$", part)
                if m and m.group(1) in code_rank:
                    names[code_rank[m.group(1)]] = m.group(2)
                    conf[code_rank[m.group(1)]] = float(m.group(3))
            parsed[p[0].split(";")[0]] = (names, conf)
    rows = []
    for asv in asv_ids(args.asvs):
        names, conf = parsed.get(asv, ({}, {}))
        keep = [conf[r] >= args.cutoff if r in conf else True for r in RANKS]
        rows.append([asv] + truncate([names.get(r, "") for r in RANKS], keep)
                    + [conf[r] if r in conf else NA for r in RANKS])
    write_table(args.out, rows, [f"conf_{r}" for r in RANKS])


def rdp(args):
    parsed = {}
    with open(args.raw) as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            names, conf, i = {}, {}, 1
            while i + 2 < len(p):
                if p[i] in ("", "-"):                       # empty or orientation column
                    i += 1
                    continue
                name, rank, c = p[i], p[i + 1], p[i + 2]
                if rank in RANKS:
                    name = re.sub(r"__\d+$", "", name)      # uniqueness suffix of the training taxonomy
                    names[rank] = "" if name.startswith("unassigned_") else name
                    conf[rank] = float(c)
                i += 3
            parsed[p[0]] = (names, conf)
    rows = []
    for asv in asv_ids(args.asvs):
        names, conf = parsed.get(asv, ({}, {}))
        keep = [conf[r] >= args.cutoff if r in conf else True for r in RANKS]
        rows.append([asv] + truncate([names.get(r, "") for r in RANKS], keep)
                    + [conf[r] if r in conf else NA for r in RANKS])
    write_table(args.out, rows, [f"conf_{r}" for r in RANKS])


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("method", choices=["blast", "sintax", "rdp"])
    ap.add_argument("--asvs", required=True)
    ap.add_argument("--raw", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--lineages")
    ap.add_argument("--identity", default="98,95,90,85,80,75")
    ap.add_argument("--cutoff", type=float, default=0.75)
    args = ap.parse_args()
    {"blast": blast, "sintax": sintax, "rdp": rdp}[args.method](args)


if __name__ == "__main__":
    main()
