#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"
FEATURE="data/COI_featureTable.qza"
TAXONOMY="data/COI_taxonomy.qza"
METADATA="metadata/COI_metadata_original.tsv"
WORK="work"; EXPORT="exports"; RESULTS="results"; CONFIDENCE="0.90"
HOLO_TERMS="Calanoida,Cyclopoida,Euphausiacea,Pteropoda,Siphonophorae,Branchiopoda,Ostracoda,Appendicularia,Thaliacea,Chaetognatha,Ctenophora,Microsetella,Euterpina,Atlanta,Hyperiidea"
mkdir -p "$WORK" "$EXPORT" "$RESULTS" "$EXPORT/beta_stats"
for f in "$FEATURE" "$TAXONOMY" "$METADATA"; do [[ -f "$f" ]] || { echo "ERROR: missing $f" >&2; exit 1; }; done
for cmd in qiime biom python; do command -v "$cmd" >/dev/null 2>&1 || { echo "ERROR: '$cmd' not found" >&2; exit 1; }; done
qiime deicode rpca --help >/dev/null 2>&1 || { echo "ERROR: DEICODE plugin unavailable" >&2; exit 1; }

qiime taxa filter-table --i-table "$FEATURE" --i-taxonomy "$TAXONOMY" --p-include "$HOLO_TERMS" --p-mode contains --o-filtered-table "$WORK/00_holoplankton_no_confidence.qza"
qiime feature-table filter-features --i-table "$FEATURE" --m-metadata-file "$TAXONOMY" --p-where "[Confidence] >= ${CONFIDENCE}" --o-filtered-table "$WORK/01_all_metazoa_conf90.qza"
qiime taxa filter-table --i-table "$WORK/01_all_metazoa_conf90.qza" --i-taxonomy "$TAXONOMY" --p-include "$HOLO_TERMS" --p-mode contains --o-filtered-table "$WORK/02_holoplankton_conf90_all20.qza"
qiime feature-table filter-samples --i-table "$WORK/02_holoplankton_conf90_all20.qza" --m-metadata-file "$METADATA" --p-where "[layer] IN ('upper', 'lower')" --o-filtered-table "$WORK/03_holoplankton_conf90_noN12.qza"

qiime deicode rpca --i-table "$WORK/03_holoplankton_conf90_noN12.qza" --p-min-feature-count 10 --p-min-sample-count 500 --o-biplot "$WORK/04_deicode_ordination.qza" --o-distance-matrix "$WORK/04_deicode_distance.qza"
for spec in "ITW:itw" "fraction:mesh" "layer:layer"; do
  col="${spec%%:*}"; tag="${spec##*:}"
  qiime diversity beta-group-significance --i-distance-matrix "$WORK/04_deicode_distance.qza" --m-metadata-file "$METADATA" --m-metadata-column "$col" --p-method permanova --p-permutations 999 --o-visualization "$EXPORT/beta_stats/permanova_${tag}.qzv"
done

rm -rf "$EXPORT/original" "$EXPORT/holo_no_conf" "$EXPORT/conf90_all20" "$EXPORT/conf90_noN12" "$EXPORT/taxonomy" "$EXPORT/ordination"
qiime tools export --input-path "$FEATURE" --output-path "$EXPORT/original"
qiime tools export --input-path "$WORK/00_holoplankton_no_confidence.qza" --output-path "$EXPORT/holo_no_conf"
qiime tools export --input-path "$WORK/02_holoplankton_conf90_all20.qza" --output-path "$EXPORT/conf90_all20"
qiime tools export --input-path "$WORK/03_holoplankton_conf90_noN12.qza" --output-path "$EXPORT/conf90_noN12"
qiime tools export --input-path "$TAXONOMY" --output-path "$EXPORT/taxonomy"
qiime tools export --input-path "$WORK/04_deicode_ordination.qza" --output-path "$EXPORT/ordination"
cp "$METADATA" "$EXPORT/COI_metadata_original.tsv"
for d in original holo_no_conf conf90_all20 conf90_noN12; do biom convert -i "$EXPORT/$d/feature-table.biom" -o "$EXPORT/$d/feature-table.tsv" --to-tsv; done

python - <<'PYCOUNT'
from biom import load_table
import csv
items=[('Original miniCOI feature table','exports/original/feature-table.biom'),('Strict holoplankton before confidence filter','exports/holo_no_conf/feature-table.biom'),('Strict holoplankton after confidence >=0.90','exports/conf90_all20/feature-table.biom'),('Statistical subset, N12 excluded','exports/conf90_noN12/feature-table.biom')]
rows=[]
for label,path in items:
    t=load_table(path); rows.append([label,t.shape[0],t.shape[1],int(t.sum())])
with open('results/filter_summary.csv','w',newline='') as f:
    w=csv.writer(f); w.writerow(['Dataset','ASVs','Samples','Reads']); w.writerows(rows)
PYCOUNT

echo "QIIME 2 filtering and beta-diversity branch complete."
