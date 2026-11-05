#!/usr/bin/env bash
# run_all.sh - taxonomic assignment of ASVs, all steps in order. Settings: config.sh.
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
source "${1:-$HERE/config.sh}"
: "${TAXDUMP_RELEASE:?set TAXDUMP_RELEASE in config.sh}"
PY="python3 $HERE/python"
mkdir -p "$OUT"/{references,raw,tables,logs}

# 0. Record of the run: software versions and settings
{ date -u +"started %Y-%m-%dT%H:%M:%SZ"
  blastn -version | sed -n 1p; vsearch --version 2>&1 | sed -n 1p; "$JAVA" -version 2>&1 | sed -n 1p; python3 --version
  sha256sum "$RDP_CLASSIFIER/classifier.jar"
  grep -E '^[A-Z_0-9]+=' "${1:-$HERE/config.sh}"; } > "$OUT/logs/run_info.txt"

# 1. NCBI taxonomy dump (dated release)
if [ ! -s "$TAXDUMP_DIR/nodes.dmp" ]; then
  mkdir -p "$TAXDUMP_DIR"
  ARCH=https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump_archive
  F=$(curl -fsS "$ARCH/" | grep -o "taxdmp_${TAXDUMP_RELEASE}[^\"<>]*\.zip" | sort -u | sed -n 1p)
  [ -n "$F" ] || { echo "No NCBI taxonomy release $TAXDUMP_RELEASE in $ARCH" >&2; exit 1; }
  curl -fsS -o "$TAXDUMP_DIR/$F" "$ARCH/$F"
  unzip -o -q "$TAXDUMP_DIR/$F" nodes.dmp names.dmp -d "$TAXDUMP_DIR"
fi

# 2. Cross-filled reference databases
$PY/build_references.py --midori "$MIDORI2" --bold "$BOLD" --mzg "$MZG" --taxdump "$TAXDUMP_DIR" \
  --out "$OUT/references" | tee "$OUT/logs/cross_fill.txt"
R="$OUT/references"

# 3. Reference indexes
makeblastdb -in "$R/blast_ref.fasta" -dbtype nucl -parse_seqids -out "$R/blast_ref" > "$OUT/logs/makeblastdb.txt"
CP="$R/java:$RDP_CLASSIFIER/classifier.jar:$RDP_CLASSIFIER/lib/*"
mkdir -p "$R/java" "$R/rdp_model"
"$JAVAC" -cp "$CP" -d "$R/java" "$HERE"/java/*.java
"$JAVA" -Xmx"$JAVA_MEM" -Dwordsize="$RDP_WORD_LENGTH" -cp "$CP" WordSizeMain train \
  -o "$R/rdp_model" -s "$R/rdp_train.fasta" -t "$R/rdp_taxonomy.txt" > "$OUT/logs/rdp_train.txt" 2>&1 || true
[ -s "$R/rdp_model/genus_wordConditionalProbList.txt" ] || { cat "$OUT/logs/rdp_train.txt" >&2; echo "RDP training failed" >&2; exit 1; }
printf '%s\n' "bergeyTree=bergeyTrainingTree.xml" "probabilityList=genus_wordConditionalProbList.txt" \
  "probabilityIndex=wordConditionalProbIndexArr.txt" "wordPrior=logWordPrior.txt" \
  "classifierVersion=MZG COI cross-filled, word length $RDP_WORD_LENGTH" > "$R/rdp_model/rRNAClassifier.properties"

# 4. Classification of each library with the three methods (a raw output that already exists is kept)
for lib in $LIBRARIES; do
  A="$ASV_DIR/$lib.ASVs.fa"
  [ -s "$OUT/raw/$lib.blast.tsv" ] || blastn -query "$A" -db "$R/blast_ref" -evalue "$BLAST_EVALUE" -max_target_seqs 1 -num_threads "$THREADS" \
    -outfmt "6 qseqid sseqid pident length qlen slen evalue bitscore" > "$OUT/raw/$lib.blast.tsv"
  [ -s "$OUT/raw/$lib.sintax.tsv" ] || vsearch --sintax "$A" --db "$R/sintax_ref.fasta" --tabbedout "$OUT/raw/$lib.sintax.tsv" \
    --sintax_cutoff "$SINTAX_CUTOFF" --randseed "$SINTAX_SEED" --threads 1 --quiet
  [ -s "$OUT/raw/$lib.rdp.tsv" ] || "$JAVA" -Xmx"$JAVA_MEM" -Dwordsize="$RDP_WORD_LENGTH" -cp "$CP" WordSizeMain classify \
    -t "$R/rdp_model/rRNAClassifier.properties" -f allrank -o "$OUT/raw/$lib.rdp.tsv" "$A"
  $PY/classification_tables.py blast  --asvs "$A" --raw "$OUT/raw/$lib.blast.tsv" --lineages "$R/ref_lineages.tsv" \
    --identity "$BLAST_IDENTITY" --out "$OUT/tables/${lib}_BLAST_MIDORI2.tsv"
  $PY/classification_tables.py sintax --asvs "$A" --raw "$OUT/raw/$lib.sintax.tsv" --cutoff "$SINTAX_CUTOFF" \
    --out "$OUT/tables/${lib}_SINTAX_BOLD.tsv"
  $PY/classification_tables.py rdp    --asvs "$A" --raw "$OUT/raw/$lib.rdp.tsv" --cutoff "$RDP_CUTOFF" \
    --out "$OUT/tables/${lib}_RDP_MZG.tsv"
done

# 5. Checksums of inputs and results
{ sha256sum "$MIDORI2" "$BOLD" "$MZG" "$TAXDUMP_DIR"/nodes.dmp "$TAXDUMP_DIR"/names.dmp
  for lib in $LIBRARIES; do sha256sum "$ASV_DIR/$lib.ASVs.fa"; done
  sha256sum "$R"/*.fasta "$R"/*.txt "$R"/*.tsv "$OUT"/raw/* "$OUT"/tables/*; } > "$OUT/logs/checksums.sha256"
date -u +"finished %Y-%m-%dT%H:%M:%SZ" >> "$OUT/logs/run_info.txt"
