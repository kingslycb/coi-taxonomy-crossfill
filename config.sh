# config.sh - settings for run_all.sh. Paths may be absolute or relative to the folder run_all.sh is started from.

# Input data: one FASTA file of ASVs per library, named <library>.ASVs.fa
ASV_DIR=asv
LIBRARIES="library1 library2 library3 library4 library5 library6 library7"

# Reference databases (downloaded files)
MIDORI2=references/MIDORI2_TOTAL_NUC_GB272_CO1_RAW.fasta.gz              # MIDORI2, GenBank release 272
BOLD=references/BOLDistilled_COI_Apr2026_SEQUENCES_sintax.fasta            # BOLDistilled COI, April 2026
MZG=references/MZGdada2-coi__MZGdbALL__o00__AC.fasta.gz                    # MetaZooGene COI, all regions, modes A + C
TAXDUMP_RELEASE=                  # NCBI taxonomy dump of the monthly archive, e.g. 2026-08-01
TAXDUMP_DIR=references/taxdump

# Software
RDP_CLASSIFIER=rdp_classifier_2.14/dist   # folder with classifier.jar and lib/
JAVA=java                         # java and javac of a JDK 21, e.g. /opt/jdk-21/bin/java
JAVAC=javac
JAVA_MEM=14g
THREADS=4

# Classification
BLAST_EVALUE=1e-50
BLAST_IDENTITY=98,95,90,85,80,75  # minimum % identity for species, genus, family, order, class, phylum
SINTAX_CUTOFF=0.75
SINTAX_SEED=1                     # fixed seed with one thread: repeatable SINTAX results
RDP_WORD_LENGTH=7
RDP_CUTOFF=0.75

# Output
OUT=results
