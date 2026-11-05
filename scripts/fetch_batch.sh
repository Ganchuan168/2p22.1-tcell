#!/bin/bash
# Batch download of GWAS Catalog immunophenotype summary statistics;
# rows with p < 1e-4 are retained and the raw file is then deleted.
# Usage: fetch_batch.sh <id_list_file> <n_parallel>
set -u
ROOT="${NPM_MR_ROOT:-$(pwd)}"
RAW=$ROOT/data/raw
OUT=$ROOT/data/ivs_filtered
mkdir -p "$RAW" "$OUT"
LIST=$1
PAR=${2:-8}

fetch_one() {
  id=$1
  num=${id#GCST}
  if [ "$num" -le 90002000 ]; then grp=GCST90001001-GCST90002000; else grp=GCST90002001-GCST90003000; fi
  url="https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/${grp}/${id}/${id}_buildGRCh37.tsv.gz"
  out="$OUT/${id}.tsv"
  [ -s "$out" ] && return 0

  tmp="$RAW/${id}.gz"
  for attempt in 1 2 3; do
    aria2c -x4 -s4 -k1M -q --allow-overwrite=true --auto-file-renaming=false \
      --console-log-level=error -d "$RAW" -o "${id}.gz" "$url" >/dev/null 2>&1
    if [ -f "$tmp" ] && [ ! -f "${tmp}.aria2" ]; then
      pigz -dc "$tmp" 2>/dev/null | awk -F'\t' 'NR==1{print "chr\tpos\tea\toa\tmaf\teaf\tbeta\tse\tp\tn"; next} $10<1e-4{print $1"\t"$2"\t"$3"\t"$4"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$5}' > "${out}.tmp"
      if [ -s "${out}.tmp" ]; then mv "${out}.tmp" "$out"; rm -f "$tmp"; return 0; fi
    fi
    rm -f "$tmp" "${tmp}.aria2"
    sleep 2
  done
  echo "FAIL $id" >> "$ROOT/data/fetch_failures.log"
  return 1
}
export -f fetch_one
export ROOT RAW OUT

cat "$LIST" | xargs -P "$PAR" -I{} bash -c 'fetch_one "$@"' _ {}
echo "done $(date)"
