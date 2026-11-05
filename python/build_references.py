#!/usr/bin/env python3
"""
Cross-fill the three COI reference databases and write them in the formats of the three classifiers.

Every record is given a seven-rank lineage (kingdom to species). Animal species (excluding arachnids and insects)
that occur in one database but not in another are copied, with their sequences and lineage, into the database that
lacks them. Species are compared by their binomial ('genus epithet').

Inputs
  --midori   MIDORI2 COI FASTA (NCBI-style headers with rank_name_taxid tokens), plain or .gz
  --bold     BOLDistilled COI FASTA in SINTAX format (;tax=k:...,s:...)
  --mzg      MetaZooGene COI FASTA in DADA2 format (semicolon-separated lineage)
  --taxdump  folder with the NCBI taxonomy dump (nodes.dmp, names.dmp); used to give MZG names their ranks
Outputs (in --out)
  blast_ref.fasta     MIDORI2 + added records, for makeblastdb
  sintax_ref.fasta    BOLDistilled + added records, for vsearch --sintax
  rdp_train.fasta     MZG + added records, with lineage, for RDP Classifier training
  rdp_taxonomy.txt    taxonomy file for RDP Classifier training
  ref_lineages.tsv    every reference record: id, database, source database, original id, seven ranks
  cross_fill_summary.tsv
"""
import argparse
import os
from collections import Counter, defaultdict

from taxutils import (RANKS, can_cross_fill, load_ncbi_ranks, mzg_tokens, clean_name, parse_midori, parse_mzg,
                      parse_sintax_header, read_fasta, sintax_label, species_key)

DATABASES = ["MIDORI2", "BOLD", "MZG"]


def records(db, args, ncbi_rank):
    """Yield (original id, lineage, sequence) for one source database."""
    if db == "MIDORI2":
        for header, seq in read_fasta(args.midori):
            rid, lineage = parse_midori(header)
            yield rid, lineage, seq
    elif db == "BOLD":
        for header, seq in read_fasta(args.bold):
            rid, lineage = parse_sintax_header(header)
            yield rid, lineage, seq
    else:
        for i, (header, seq) in enumerate(read_fasta(args.mzg), 1):
            yield f"MZG{i}", parse_mzg(header, ncbi_rank), seq


class RdpTree:
    """Lineages -> RDP taxonomy file. The RDP Classifier requires names to be unique within a rank, compared
    without regard to case, so a name that occurs under different parents gets a suffix '__2', '__3', ...
    (removed again when the classifier output is read)."""

    def __init__(self):
        self.nodes = {(): (0, "Root", -1, 0, "rootrank")}          # path -> (taxid, name, parent taxid, depth, rank)
        self.used = defaultdict(set)                               # rank -> names in use, lower case

    def lineage_names(self, lineage):
        """Seven names for the training header; missing ranks become placeholders; no genus -> None."""
        if "genus" not in lineage:
            return None
        names, path = [], ()
        for depth, rank in enumerate(RANKS, 1):
            if rank == "species":
                name = lineage.get("species") or f"{lineage['genus']} sp."
            else:
                name = lineage.get(rank) or f"unassigned_{rank}"
            key = path + (name,)
            if key not in self.nodes:
                unique, k = name, 1
                while unique.lower() in self.used[rank]:
                    k += 1
                    unique = f"{name}__{k}"
                self.used[rank].add(unique.lower())
                self.nodes[key] = (len(self.nodes), unique, self.nodes[path][0], depth, rank)
            names.append(self.nodes[key][1])
            path = key
        return names

    def write(self, path):
        with open(path, "w") as f:
            for taxid, name, parent, depth, rank in sorted(self.nodes.values()):
                f.write(f"{taxid}*{name}*{parent}*{depth}*{rank}\n")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    for a in ("midori", "bold", "mzg", "taxdump", "out"):
        ap.add_argument(f"--{a}", required=True)
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    # NCBI ranks, only for the names that occur in MZG lineages
    wanted = set()
    for header, _ in read_fasta(args.mzg):
        wanted.update(clean_name(t) for t in mzg_tokens(header))
    ncbi_rank = load_ncbi_ranks(os.path.join(args.taxdump, "nodes.dmp"), os.path.join(args.taxdump, "names.dmp"),
                                wanted)

    # pass 1: species that can be cross-filled, per database
    species = {db: set() for db in DATABASES}
    for db in DATABASES:
        for _, lineage, _ in records(db, args, ncbi_rank):
            if can_cross_fill(lineage):
                species[db].add(species_key(lineage))
    union = set().union(*species.values())
    missing = {db: union - species[db] for db in DATABASES}

    # pass 2: each database = its own records + records of its missing species from the other two
    out = {"MIDORI2": open(os.path.join(args.out, "blast_ref.fasta"), "w"),
           "BOLD": open(os.path.join(args.out, "sintax_ref.fasta"), "w"),
           "MZG": open(os.path.join(args.out, "rdp_train.fasta"), "w")}
    lineages = open(os.path.join(args.out, "ref_lineages.tsv"), "w")
    lineages.write("ref_id\tdatabase\tsource_database\toriginal_id\t" + "\t".join(RANKS) + "\n")
    tree, n, counts, skipped_rdp = RdpTree(), 0, Counter(), 0

    def write(target, source, rid, lineage, seq):
        nonlocal n, skipped_rdp
        if target == "MZG":
            names = tree.lineage_names(lineage)
            if names is None:
                skipped_rdp += 1
                return
        n += 1
        ref_id = f"r{n:08d}"
        if target == "MIDORI2":
            out[target].write(f">{ref_id}\n{seq}\n")
        elif target == "BOLD":
            out[target].write(f">{ref_id};tax={sintax_label(lineage)};\n{seq}\n")
        else:
            out[target].write(f">{ref_id}\tRoot;{';'.join(names)}\n{seq}\n")
        lineages.write("\t".join([ref_id, target, source, rid] + [lineage.get(r, "") for r in RANKS]) + "\n")
        counts[(target, source)] += 1

    for source in DATABASES:
        for rid, lineage, seq in records(source, args, ncbi_rank):
            write(source, source, rid, lineage, seq)
            if can_cross_fill(lineage):
                key = species_key(lineage)
                for target in DATABASES:
                    if target != source and key in missing[target]:
                        write(target, source, rid, lineage, seq)
    for f in list(out.values()) + [lineages]:
        f.close()
    tree.write(os.path.join(args.out, "rdp_taxonomy.txt"))

    lines = ["database\tspecies_before\tspecies_added\trecords_own\trecords_added"]
    for db in DATABASES:
        added = sum(v for (t, s), v in counts.items() if t == db and s != db)
        lines.append(f"{db}\t{len(species[db])}\t{len(missing[db])}\t{counts[(db, db)]}\t{added}")
    with open(os.path.join(args.out, "cross_fill_summary.tsv"), "w") as f:
        f.write("\n".join(lines) + "\n")
    print("\n".join(lines))
    if skipped_rdp:
        print(f"records without a genus, not used for RDP training: {skipped_rdp}")


if __name__ == "__main__":
    main()
