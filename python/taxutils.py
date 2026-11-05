"""Shared helpers: FASTA reading, header parsing for the three reference databases, seven-rank lineages."""
import gzip
import re

RANKS = ["kingdom", "phylum", "class", "order", "family", "genus", "species"]
SINTAX_CODE = dict(zip(RANKS, "kpcofgs"))
ANIMAL_KINGDOMS = {"Animalia", "Metazoa"}
EXCLUDED_CLASSES = {"Arachnida", "Insecta"}          # arachnids and insects are not cross-filled
NOT_EPITHETS = {"sp", "spp", "cf", "aff", "nr", "gen", "indet"}


def open_text(path):
    """Open a plain or gzip-compressed text file for reading."""
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def read_fasta(path):
    """Yield (header without '>', sequence) for every record; handles wrapped sequences."""
    header, seq = None, []
    with open_text(path) as f:
        for line in f:
            line = line.rstrip("\n\r")
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq)
                header, seq = line[1:], []
            elif line:
                seq.append(line.strip())
    if header is not None:
        yield header, "".join(seq)


def clean_name(name):
    """Name as used in all outputs: underscores to spaces, separators used by the file formats removed."""
    name = re.sub(r"[_\s]+", " ", name).strip()
    return re.sub(r"[,;:*\t]", "", name)


def drop_subgenus(name):
    """'Acartia (Acartiura) clausi' or 'Acartia Acartiura clausi' -> 'Acartia clausi'. The subgenus is not part of
    the binomial and is written differently (or not at all) in the three databases. Qualifiers such as 'Cf.' or
    'Aff.' are kept, so an uncertain identification stays uncertain."""
    name = re.sub(r"\s*\([^)]*\)\s*", " ", name).strip()
    words = name.split()
    if len(words) >= 3 and words[1][:1].isupper() and words[1].lower().rstrip(".") not in NOT_EPITHETS:
        words = [words[0]] + words[2:]
    return " ".join(words)


def normalise(lineage):
    """Seven-rank lineage with one label for the animal kingdom and species names without subgenus;
    empty ranks are dropped."""
    out = {r: clean_name(lineage[r]) for r in RANKS if lineage.get(r)}
    if "species" in out:
        out["species"] = drop_subgenus(out["species"])
    out = {r: v for r, v in out.items() if v}
    if out.get("kingdom") in ANIMAL_KINGDOMS:
        out["kingdom"] = "Animalia"
    if "species" in out and "genus" not in out:
        out["genus"] = out["species"].split()[0]
    return out


def species_key(lineage):
    """'genus epithet' in lower case for a species-level name, or None (sp., cf., aff. and similar are not species)."""
    words = lineage.get("species", "").split()
    if len(words) < 2:
        return None
    genus, epithet = words[0], words[1].rstrip(".")
    if not re.fullmatch(r"[A-Z][a-z]+", genus) or not re.fullmatch(r"[a-z][a-z-]+", epithet) or epithet in NOT_EPITHETS:
        return None
    return f"{genus.lower()} {epithet}"


def can_cross_fill(lineage):
    """Animal species-level records, excluding arachnids and insects."""
    return (lineage.get("kingdom") == "Animalia" and lineage.get("class") not in EXCLUDED_CLASSES
            and species_key(lineage) is not None)


# --- MIDORI2: "<accession> root_1;...;kingdom_Metazoa_33208;phylum_Arthropoda_6656;...;species_Calanus_finmarchicus_6837"
_MIDORI_TOKEN = re.compile(r"^(kingdom|phylum|class|order|family|genus|species)_(.+?)(?:_\d+)?$")


def parse_midori(header):
    """Names in the rank tokens may contain spaces or underscores; the accession is the first word."""
    parts = header.split(None, 1)
    accession = parts[0]
    lineage = {}
    for token in (parts[1] if len(parts) > 1 else "").split(";"):
        m = _MIDORI_TOKEN.match(token.strip())
        if m:
            lineage[m.group(1)] = m.group(2)
    return accession, normalise(lineage)


# --- BOLDistilled (SINTAX format): "<id>[|...];tax=k:Animalia,p:Arthropoda,...,s:Calanus finmarchicus;"
def parse_sintax_header(header):
    record_id = re.split(r"[;|\s]", header, 1)[0]
    m = re.search(r"tax=([^;]*)", header)
    lineage = {}
    if m:
        codes = {v: k for k, v in SINTAX_CODE.items()}
        for part in m.group(1).split(","):
            code, _, name = part.partition(":")
            if code in codes and name:
                lineage[codes[code]] = name
    return record_id, normalise(lineage)


# --- MetaZooGene (DADA2 format): "Animalia;Animalia_EXT;Arthropoda;Crustacea;...;Calanidae;Calanus;Calanus_finmarchicus;"
def mzg_tokens(header):
    """Lineage names without the '_EXT' placeholder levels."""
    return [t for t in header.strip().rstrip(";").split(";") if t and not t.endswith("_EXT")]


def parse_mzg(header, ncbi_rank):
    """Seven-rank lineage of an MZG record. The kingdom is the first name; phylum to family come from the rank
    that the NCBI taxonomy gives each name (names with no or several canonical ranks are skipped; ranks must follow
    the lineage order); genus and species come from the last name, which MZG gives as 'Genus_species'."""
    tokens = mzg_tokens(header)
    if not tokens:
        return {}
    lineage = {"kingdom": tokens[0]}
    last = clean_name(tokens[-1])
    if len(last.split()) >= 2:
        lineage["species"] = last
        lineage["genus"] = last.split()[0]
        middle = tokens[1:-1]
    else:
        middle = tokens[1:]
    level = 0
    for token in middle:
        ranks = ncbi_rank.get(clean_name(token), set()) & set(RANKS[1:6])
        if len(ranks) == 1:
            rank = ranks.pop()
            i = RANKS.index(rank)
            if i > level and rank not in lineage:
                lineage[rank] = token
                level = i
    return normalise(lineage)


def load_ncbi_ranks(nodes_dmp, names_dmp, wanted):
    """name -> set of NCBI ranks, for the scientific names in `wanted` only."""
    rank_of = {}
    with open(nodes_dmp) as f:
        for line in f:
            p = line.split("\t|\t")
            rank_of[p[0]] = p[2]
    out = {}
    with open(names_dmp) as f:
        for line in f:
            p = [x.strip() for x in line.rstrip("\t|\n").split("\t|\t")]
            if len(p) >= 4 and p[3] == "scientific name" and p[1] in wanted:
                out.setdefault(p[1], set()).add(rank_of.get(p[0], ""))
    return out


def sintax_label(lineage):
    return ",".join(f"{SINTAX_CODE[r]}:{lineage[r]}" for r in RANKS if r in lineage)
