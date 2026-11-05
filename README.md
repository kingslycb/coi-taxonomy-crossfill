# COI metabarcoding: taxonomic assignment of ASVs

Taxonomic assignment of amplicon sequence variants (ASVs) against three cross-filled COI reference databases with
three classifiers. ASVs were generated beforehand with Cutadapt and DADA2 (see the article).

## Workflow

1. **Reference databases** (`python/build_references.py`).
   - Inputs: MIDORI2 (GenBank release 272), BOLDistilled (April 2026) and MetaZooGene (MZGdb COI, all regions,
     modes A and C).
   - Every record gets a seven-rank lineage (kingdom to species). MZG lineages have no rank labels; their ranks are
     taken from the NCBI taxonomy (a dated monthly release).
   - Species names are reduced to the binomial: a subgenus, written with or without parentheses
     ("Acartia (Acartiura) clausi", "Acartia Acartiura clausi"), is removed, because the databases write it
     differently or not at all. Qualifiers (sp., cf., aff.) are kept, and such names are not treated as species.
   - Animal species (excluding arachnids and insects) found in one database but not in another are copied into the
     database that lacks them, with their sequences and lineages.
2. **Classification**, each ASV separately against each database (`run_all.sh`):
   - BLASTn against MIDORI2: E-value ≤ 1e-50, top hit. The hit's lineage is kept down to the finest rank whose
     identity threshold is met: 98% species, 95% genus, 90% family, 85% order, 80% class, 75% phylum.
   - SINTAX (VSEARCH) against BOLDistilled: bootstrap cut-off 0.75.
   - RDP naive Bayesian classifier trained on MZG, word length 7 (`java/WordSizeMain.java`; `java/TreeFactory.java`
     is the RDP training class with its word arrays sized for word lengths other than 8). Ranks with
     confidence below 0.75 are left unassigned. Reference sequences without a species name are trained as
     "Genus sp.".
3. **Result tables** (`python/classification_tables.py`). One table per library and method, with the same
   columns: `ASV`, `kingdom` … `species` (`NA` = unassigned), plus the method's support values. Once a rank fails
   its threshold, all finer ranks are unassigned.

## Requirements

- The software in `environment.yml`:
  - BLAST+ 2.17.0
  - VSEARCH 2.30.1
  - Java (JDK) 21
  - Python 3
- RDP Classifier 2.14 (SourceForge project `rdp-classifier`). `java/TreeFactory.java` is a modified copy of one of its classes
  (GNU General Public License, version 2 or any later version); the changes are listed at the top of the file.

```bash
conda env create -f environment.yml
conda activate coi-taxonomy
```

## Input data

- ASVs: one FASTA file per sequencing library (`<library>.ASVs.fa`)
- `MIDORI2_TOTAL_NUC_GB272_CO1_RAW.fasta.gz`: MIDORI2 (Leray et al., 2022)
- `BOLDistilled_COI_Apr2026_SEQUENCES_sintax.fasta`: BOLDistilled (Prosser et al., 2025)
- `MZGdada2-coi__MZGdbALL__o00__AC.fasta.gz`: MetaZooGene (Bucklin et al., 2021)
- The NCBI taxonomy release is downloaded by `run_all.sh` from
  https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump_archive/

## Running

1. Set the paths, library names and NCBI taxonomy release in `config.sh`.
2. Run:

```bash
bash run_all.sh            # or: bash run_all.sh my_config.sh
```

The reference databases are always rebuilt. Classifier outputs that already exist in `results/raw/` are kept, so
a run that stopped can be restarted; delete `results/` for a complete new run.

## Outputs (`results/`)

| Folder | Content |
|---|---|
| `references/` | cross-filled databases, `ref_lineages.tsv` (every reference record and its source), `cross_fill_summary.tsv`, BLAST and RDP indexes |
| `raw/` | raw output of BLASTn, SINTAX and the RDP Classifier |
| `tables/` | `<library>_BLAST_MIDORI2.tsv`, `<library>_SINTAX_BOLD.tsv`, `<library>_RDP_MZG.tsv` |
| `logs/` | software versions and settings of the run (`run_info.txt`), SHA-256 checksums of inputs and results |

## Reproducibility

- **Software:** versions are fixed in `environment.yml` and recorded in `logs/run_info.txt` at each run.
- **Taxonomy:** the NCBI taxonomy is a dated release, set in `config.sh`.
- **Random numbers:**
  - SINTAX draws random subsamples. It runs with a fixed seed and one thread, which gives the same result every
    time.
  - The RDP Classifier chooses words with a fixed seed. It breaks ties between equally supported taxa at random.
- **Ambiguous names:** reference record IDs and RDP taxon names are assigned in input order, so the same inputs
  always give the same files.

## Tests

```bash
python3 tests/test_pipeline.py
```

The tests use small made-up files. They cover cross-filling, the seven-rank lineages, the RDP training files and
the three result tables, and they run SINTAX when VSEARCH is installed.

## Licence and citation

GNU General Public License, version 3 or later (`LICENSE`). `java/TreeFactory.java` is derived from the RDP
Classifier (GNU General Public License, version 2 or any later version). To cite this software, see `CITATION.cff`.

## References

- Bucklin, A., Peijnenburg, K. T. C. A., Kosobokova, K. N., O'Brien, T. D., Blanco-Bercial, L., Cornils, A., Falkenhaug, T., Hopcroft, R. R., Hosia, A., Laakmann, S., Li, C., Martell, L., Questel, J. M., Wall-Palmer, D., Wang, M., Wiebe, P. H., & Weydmann-Zwolicka, A. (2021). Toward a global reference database of COI barcodes for marine zooplankton. *Marine Biology*, 168(6), 78. https://doi.org/10.1007/s00227-021-03887-y
- Camacho, C., Coulouris, G., Avagyan, V., Ma, N., Papadopoulos, J., Bealer, K., & Madden, T. L. (2009). BLAST+: architecture and applications. *BMC Bioinformatics*, 10, 421. https://doi.org/10.1186/1471-2105-10-421
- Edgar, R. C. (2016). SINTAX: a simple non-Bayesian taxonomy classifier for 16S and ITS sequences. *bioRxiv*, 074161. https://doi.org/10.1101/074161
- Leray, M., Knowlton, N., & Machida, R. J. (2022). MIDORI2: A collection of quality controlled, preformatted, and regularly updated reference databases for taxonomic assignment of eukaryotic mitochondrial sequences. *Environmental DNA*, 4(4), 894–907. https://doi.org/10.1002/edn3.303
- Prosser, S. W. J., Floyd, R. M., Thompson, K. A., & Hebert, P. D. N. (2025). BOLDistilled: Comprehensive but compact DNA barcode reference libraries. *Molecular Ecology Resources*. https://doi.org/10.1111/1755-0998.70043
- Rognes, T., Flouri, T., Nichols, B., Quince, C., & Mahé, F. (2016). VSEARCH: a versatile open source tool for metagenomics. *PeerJ*, 4, e2584. https://doi.org/10.7717/peerj.2584
- Wang, Q., Garrity, G. M., Tiedje, J. M., & Cole, J. R. (2007). Naive Bayesian classifier for rapid assignment of rRNA sequences into the new bacterial taxonomy. *Applied and Environmental Microbiology*, 73(16), 5261–5267. https://doi.org/10.1128/AEM.00062-07
