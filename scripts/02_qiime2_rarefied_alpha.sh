#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"
METADATA="metadata/COI_metadata_original.tsv"; WORK="work"; EXPORT="exports"; RESULTS="results"; DEPTH=19000
TABLE="$WORK/03_holoplankton_conf90_noN12.qza"
[[ -f "$TABLE" ]] || { echo "ERROR: $TABLE not found. Run scripts/01_qiime2_filter_beta.sh first." >&2; exit 1; }
for cmd in qiime python; do command -v "$cmd" >/dev/null 2>&1 || { echo "ERROR: '$cmd' not found" >&2; exit 1; }; done
mkdir -p "$EXPORT/alpha_table" "$EXPORT/alpha_observed" "$EXPORT/alpha_shannon" "$RESULTS"
rm -rf "$EXPORT/alpha_table"/* "$EXPORT/alpha_observed"/* "$EXPORT/alpha_shannon"/*
qiime tools export --input-path "$TABLE" --output-path "$EXPORT/alpha_table"
MIN_READS="$(python - "$EXPORT/alpha_table/feature-table.biom" <<'PYDEPTH'
import sys
from biom import load_table
t=load_table(sys.argv[1]); print(int(min(t.sum(axis='sample'))))
PYDEPTH
)"
if (( MIN_READS < DEPTH )); then echo "ERROR: minimum retained depth is $MIN_READS, below manuscript depth $DEPTH." >&2; exit 1; fi
printf 'Minimum retained reads/sample before rarefaction: %s\nRarefaction depth used: %s\n' "$MIN_READS" "$DEPTH" > "$RESULTS/rarefaction_depth.txt"
qiime diversity alpha-rarefaction --i-table "$TABLE" --p-max-depth "$DEPTH" --m-metadata-file "$METADATA" --o-visualization "$RESULTS/alpha_rarefaction.qzv"
qiime feature-table rarefy --i-table "$TABLE" --p-sampling-depth "$DEPTH" --o-rarefied-table "$WORK/05_holoplankton_conf90_noN12_rarefied_19000.qza"
qiime diversity alpha --i-table "$WORK/05_holoplankton_conf90_noN12_rarefied_19000.qza" --p-metric observed_features --o-alpha-diversity "$WORK/06_observed_features_rarefied.qza"
qiime diversity alpha --i-table "$WORK/05_holoplankton_conf90_noN12_rarefied_19000.qza" --p-metric shannon --o-alpha-diversity "$WORK/06_shannon_rarefied.qza"
qiime tools export --input-path "$WORK/06_observed_features_rarefied.qza" --output-path "$EXPORT/alpha_observed"
qiime tools export --input-path "$WORK/06_shannon_rarefied.qza" --output-path "$EXPORT/alpha_shannon"
echo "Rarefied alpha-diversity branch complete at ${DEPTH} reads/sample (minimum before rarefaction: ${MIN_READS})."
