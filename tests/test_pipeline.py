#!/usr/bin/env python3
"""
Tests on small made-up files: cross-filling, seven-rank lineages, RDP training files and the three result tables.
Run: python3 tests/test_pipeline.py      (the SINTAX test also runs vsearch if it is on the PATH)
"""
import os
import random
import shutil
import subprocess
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
PY = os.path.join(HERE, "..", "python")
sys.path.insert(0, PY)
import taxutils  # noqa: E402

random.seed(1)
SEQ = {k: "".join(random.choice("ACGT") for _ in range(313)) for k in
       ["ostrea", "calanus", "insect", "diatom", "acartia", "temora", "acartia_bold"]}


def mutate(seq, n):
    s = list(seq)
    for i in random.sample(range(len(s)), n):
        s[i] = {"A": "C", "C": "G", "G": "T", "T": "A"}[s[i]]
    return "".join(s)


def write(path, text):
    with open(path, "w") as f:
        f.write(text)


def make_inputs(d):
    mid = ("root_1;cellularroot_cellular_organisms_131567;domain_Eukaryota_2759;kingdom_Metazoa_33208;"
           "phylum_{p};class_{c};order_{o};family_{f};genus_{g};species_{s}")
    write(f"{d}/midori.fasta",
          f">MT593041.1.<1.>456 " + mid.format(p="Mollusca_6447", c="Bivalvia_6544", o="Ostreida_1", f="Ostreidae_6563",
                                                g="Ostrea_6564", s="Ostrea_edulis_37623") + f"\n{SEQ['ostrea']}\n"
          f">AB1.1 " + mid.format(p="Arthropoda_6656", c="Hexanauplia_72037", o="Calanoida_6833", f="Calanidae_6834",
                                  g="Calanus_6835", s="Calanus_finmarchicus_6837") + f"\n{SEQ['calanus']}\n"
          f">AB2.1 " + mid.format(p="Arthropoda_6656", c="Insecta_50557", o="Diptera_7147", f="Culicidae_7157",
                                  g="Aedes_7158", s="Aedes_aegypti_7159") + f"\n{SEQ['insect']}\n"
          f">AB3.1 root_1;domain_Eukaryota_2759;phylum_Bacillariophyta_2836;class_Coscinodiscophyceae_2;"
          f"order_Rhizosoleniales_3;family_Rhizosoleniaceae_4;genus_Proboscia_5;species_Proboscia_alata_6\n"
          f"{SEQ['diatom']}\n")
    write(f"{d}/bold.fasta",
          f">BOLD1|X;tax=k:Animalia,p:Arthropoda,c:Hexanauplia,o:Calanoida,f:Calanidae,g:Calanus,s:Calanus finmarchicus;\n"
          f"{SEQ['calanus']}\n"
          f">BOLD2|X;tax=k:Animalia,p:Arthropoda,c:Hexanauplia,o:Calanoida,f:Acartiidae,g:Acartia,s:Acartia clausi;\n"
          f"{SEQ['acartia_bold']}\n"
          f">BOLD3|X;tax=k:Animalia,p:Arthropoda,c:Insecta,o:Diptera,f:Chironomidae,g:Chironomus,s:Chironomus riparius;\n"
          f"{SEQ['insect']}\n")
    write(f"{d}/mzg.fasta",
          f">Animalia;Animalia_EXT;Arthropoda;Crustacea;Hexanauplia;Calanoida;Acartiidae;Acartia;Acartia_clausi;\n"
          f"{SEQ['acartia'][:150]}\n{SEQ['acartia'][150:]}\n"
          f">Animalia;Arthropoda;Crustacea;Hexanauplia;Calanoida;Temoridae;Temora;Temora_longicornis;\n{SEQ['temora']}\n"
          f">Chromista;Chromista_EXT;Haptophyta;Prymnesiophyceae;Phaeocystales;Phaeocystaceae;Phaeocystis;"
          f"Phaeocystis_globosa;\n{SEQ['diatom']}\n")
    os.makedirs(f"{d}/taxdump")
    nodes = [("6656", "phylum"), ("6657", "subphylum"), ("72037", "class"), ("6833", "order"), ("6838", "family"),
             ("6839", "family"), ("2830", "phylum"), ("2831", "class"), ("2832", "order"), ("2833", "family")]
    names = [("6656", "Arthropoda"), ("6657", "Crustacea"), ("72037", "Hexanauplia"), ("6833", "Calanoida"),
             ("6838", "Acartiidae"), ("6839", "Temoridae"), ("2830", "Haptophyta"), ("2831", "Prymnesiophyceae"),
             ("2832", "Phaeocystales"), ("2833", "Phaeocystaceae")]
    write(f"{d}/taxdump/nodes.dmp", "".join(f"{t}\t|\t1\t|\t{r}\t|\t\n" for t, r in nodes))
    write(f"{d}/taxdump/names.dmp", "".join(f"{t}\t|\t{n}\t|\t\t|\tscientific name\t|\n" for t, n in names))


def read_tsv(path):
    with open(path) as f:
        cols = f.readline().rstrip("\n").split("\t")
        return [dict(zip(cols, line.rstrip("\n").split("\t"))) for line in f]


class Pipeline(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.d = tempfile.mkdtemp()
        make_inputs(cls.d)
        subprocess.run([sys.executable, f"{PY}/build_references.py", "--midori", f"{cls.d}/midori.fasta",
                        "--bold", f"{cls.d}/bold.fasta", "--mzg", f"{cls.d}/mzg.fasta", "--taxdump",
                        f"{cls.d}/taxdump", "--out", f"{cls.d}/ref"], check=True, capture_output=True)
        cls.lin = read_tsv(f"{cls.d}/ref/ref_lineages.tsv")

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.d)

    def species_in(self, db):
        return {r["species"] for r in self.lin if r["database"] == db}

    def test_parsers(self):
        self.assertEqual(taxutils.species_key({"species": "Calanus sp. A"}), None)
        self.assertEqual(taxutils.species_key({"species": "Acartia clausi"}), "acartia clausi")
        for name in ("Acartia Acartiura clausi", "Acartia (Acartiura) clausi", "Acartia_Acartiura_clausi"):
            self.assertEqual(taxutils.normalise({"species": name})["species"], "Acartia clausi")
        self.assertEqual(taxutils.normalise({"species": "Calanus sp. A"})["species"], "Calanus sp. A")
        self.assertEqual(taxutils.normalise({"species": "Calanus Cf. finmarchicus"})["species"], "Calanus Cf. finmarchicus")
        self.assertIsNone(taxutils.species_key(taxutils.normalise({"species": "Calanus Cf. finmarchicus"})))
        self.assertEqual(taxutils.parse_midori("X root_1;kingdom_Metazoa_33208;species_Ostrea_edulis_37623")[1]["kingdom"],
                         "Animalia")
        for h in ("X root_1;kingdom_Metazoa_33208;genus_Calanus_6835;species_Calanus finmarchicus_6837",
                  "X root_1;kingdom_Metazoa_33208;genus_Calanus_6835;species_Calanus_finmarchicus_6837"):
            lin = taxutils.parse_midori(h)[1]
            self.assertEqual(lin["species"], "Calanus finmarchicus")
            self.assertEqual(taxutils.species_key(lin), "calanus finmarchicus")

    def test_cross_fill(self):
        self.assertIn("Ostrea edulis", self.species_in("BOLD"))         # MIDORI2 -> BOLD
        self.assertIn("Ostrea edulis", self.species_in("MZG"))          # MIDORI2 -> MZG
        self.assertIn("Temora longicornis", self.species_in("MIDORI2"))  # MZG -> MIDORI2
        self.assertIn("Acartia clausi", self.species_in("MIDORI2"))      # BOLD and MZG -> MIDORI2
        self.assertIn("Calanus finmarchicus", self.species_in("MZG"))
        for db in ("BOLD", "MZG"):                                        # insects are not copied
            self.assertNotIn("Aedes aegypti", self.species_in(db))
        for db in ("MIDORI2", "MZG"):
            self.assertNotIn("Chironomus riparius", self.species_in(db))
        self.assertNotIn("Proboscia alata", self.species_in("BOLD"))     # non-animals are not copied

    def test_seven_ranks(self):
        for r in self.lin:
            if r["source_database"] == "MZG":
                self.assertIn(r["kingdom"], {"Animalia", "Chromista"})
                for rank in ("phylum", "class", "order", "family", "genus", "species"):
                    self.assertTrue(r[rank], (r["species"], rank))
            if r["kingdom"]:
                self.assertNotEqual(r["kingdom"], "Metazoa")
        temora = [r for r in self.lin if r["species"] == "Temora longicornis" and r["database"] == "BOLD"][0]
        self.assertEqual((temora["phylum"], temora["class"], temora["family"]), ("Arthropoda", "Hexanauplia", "Temoridae"))

    def test_mzg_non_animals(self):
        phaeo = [r for r in self.lin if r["species"] == "Phaeocystis globosa"]
        self.assertEqual({r["database"] for r in phaeo}, {"MZG"})                 # used for RDP, not cross-filled
        self.assertEqual((phaeo[0]["kingdom"], phaeo[0]["phylum"], phaeo[0]["family"]),
                         ("Chromista", "Haptophyta", "Phaeocystaceae"))

    def test_rdp_files(self):
        with open(f"{self.d}/ref/rdp_taxonomy.txt") as f:
            tax = [l.rstrip("\n").split("*") for l in f]
        ids = {t[0] for t in tax}
        self.assertEqual(tax[0], ["0", "Root", "-1", "0", "rootrank"])
        for t in tax[1:]:
            self.assertIn(t[2], ids)                                       # every parent exists
        names = [(t[4], t[1]) for t in tax]
        self.assertEqual(len(names), len(set(names)))                      # names unique within a rank
        with open(f"{self.d}/ref/rdp_train.fasta") as f:
            headers = [l for l in f if l.startswith(">")]
        for h in headers:
            self.assertEqual(len(h.split("\t")[1].strip().split(";")), 8)  # Root + seven ranks

    def test_rdp_names_unique_ignoring_case(self):
        sys.path.insert(0, PY)
        from build_references import RdpTree
        tree = RdpTree()
        a = tree.lineage_names({"kingdom": "Animalia", "phylum": "Nemertea", "genus": "Ototyphlonemertes",
                                "species": "Ototyphlonemertes duplex"})
        b = tree.lineage_names({"kingdom": "Animalia", "phylum": "Annelida", "genus": "ototyphlonemertes",
                                "species": "ototyphlonemertes duplex"})
        self.assertEqual(a[-1], "Ototyphlonemertes duplex")
        self.assertEqual(b[-1], "ototyphlonemertes duplex__2")
        names = [(n[4], n[1].lower()) for n in tree.nodes.values()]
        self.assertEqual(len(names), len(set(names)))

    def test_blast_table(self):
        lin = {r["ref_id"]: r for r in self.lin if r["database"] == "MIDORI2"}
        ostrea = [k for k, v in lin.items() if v["species"] == "Ostrea edulis"][0]
        calanus = [k for k, v in lin.items() if v["species"] == "Calanus finmarchicus"][0]
        write(f"{self.d}/q.fa", ">ASV_1\nACGT\n>ASV_2\nACGT\n>ASV_3\nACGT\n")
        write(f"{self.d}/q.blast.tsv", f"ASV_1\t{ostrea}\t99.5\t313\t313\t313\t1e-150\t560\n"
                                        f"ASV_2\t{calanus}\t91.0\t313\t313\t313\t1e-120\t420\n")
        subprocess.run([sys.executable, f"{PY}/classification_tables.py", "blast", "--asvs", f"{self.d}/q.fa",
                        "--raw", f"{self.d}/q.blast.tsv", "--lineages", f"{self.d}/ref/ref_lineages.tsv",
                        "--out", f"{self.d}/q_blast.tsv"], check=True)
        t = {r["ASV"]: r for r in read_tsv(f"{self.d}/q_blast.tsv")}
        self.assertEqual(t["ASV_1"]["species"], "Ostrea edulis")
        self.assertEqual((t["ASV_2"]["family"], t["ASV_2"]["genus"], t["ASV_2"]["species"]), ("Calanidae", "NA", "NA"))
        self.assertEqual(t["ASV_3"]["phylum"], "NA")

    def test_rdp_table(self):
        write(f"{self.d}/r.fa", ">ASV_1\nACGT\n>ASV_2\nACGT\n")
        line = ("ASV_1\t\tRoot\trootrank\t1.0\tAnimalia\tkingdom\t1.0\tArthropoda\tphylum\t1.0\tHexanauplia\tclass\t0.99\t"
                "Calanoida\torder\t0.98\tunassigned_family__2\tfamily\t0.97\tAcartia__2\tgenus\t0.9\t"
                "Acartia clausi\tspecies\t0.6\n")
        write(f"{self.d}/r.rdp.tsv", line)
        subprocess.run([sys.executable, f"{PY}/classification_tables.py", "rdp", "--asvs", f"{self.d}/r.fa",
                        "--raw", f"{self.d}/r.rdp.tsv", "--out", f"{self.d}/r_rdp.tsv"], check=True)
        t = {r["ASV"]: r for r in read_tsv(f"{self.d}/r_rdp.tsv")}
        self.assertEqual((t["ASV_1"]["family"], t["ASV_1"]["genus"], t["ASV_1"]["species"]), ("NA", "Acartia", "NA"))
        self.assertEqual(t["ASV_2"]["kingdom"], "NA")

    @unittest.skipUnless(shutil.which("vsearch"), "vsearch not on PATH")
    def test_sintax(self):
        write(f"{self.d}/s.fa", f">ASV_1\n{mutate(SEQ['ostrea'], 3)}\n>ASV_2\n{mutate(SEQ['temora'], 2)}\n")
        subprocess.run(["vsearch", "--sintax", f"{self.d}/s.fa", "--db", f"{self.d}/ref/sintax_ref.fasta",
                        "--tabbedout", f"{self.d}/s.sintax.tsv", "--sintax_cutoff", "0.75", "--randseed", "1",
                        "--threads", "1", "--quiet"], check=True)
        subprocess.run([sys.executable, f"{PY}/classification_tables.py", "sintax", "--asvs", f"{self.d}/s.fa",
                        "--raw", f"{self.d}/s.sintax.tsv", "--out", f"{self.d}/s_sintax.tsv"], check=True)
        t = {r["ASV"]: r for r in read_tsv(f"{self.d}/s_sintax.tsv")}
        self.assertEqual(t["ASV_1"]["species"], "Ostrea edulis")       # found only because it was cross-filled
        self.assertEqual(t["ASV_2"]["family"], "Temoridae")             # MZG record with all seven ranks


if __name__ == "__main__":
    unittest.main(verbosity=2)
