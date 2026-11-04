#!/bin/bash
# Build the linkage-disequilibrium matrix used by the variant-level and
# correlation-aware analyses (scripts/variant_table.R, scripts/correlated_ivw.R).
#
# Input : tables/Supplementary_Table_S1_instrument_level.csv  (produced by make_supp_table.R)
# Output: data/ld/all759.bim        variant list in PLINK .bim order
#         data/ld/all759ld.ld       753 x 753 r (not r2) matrix, tab-delimited, no header
#
# Requires PLINK 1.9 and a binary reference panel of 1000 Genomes phase 3
# European samples. Point PLINK_1KG_EUR at the prefix of that panel
# (the three files <prefix>.bed/.bim/.fam must exist).
#
# Usage: NPM_MR_ROOT=. PLINK_1KG_EUR=/path/to/EUR bash scripts/make_ld_matrix.sh
set -eu

ROOT="${NPM_MR_ROOT:-$(pwd)}"
PLINK="${PLINK_BIN:-plink1.9}"
BFILE="${PLINK_1KG_EUR:-$ROOT/data/ld/EUR}"
mkdir -p "$ROOT/data/ld"

python3 - "$ROOT" <<'PY'
import csv, os, sys
root = sys.argv[1]
ids = set()
with open(os.path.join(root, "tables/Supplementary_Table_S1_instrument_level.csv")) as fh:
    for row in csv.DictReader(fh):
        ids.add(row["rsid"])
ids = sorted(ids)
with open(os.path.join(root, "data/ld/all759.snplist"), "w") as fh:
    fh.write("\n".join(ids) + "\n")
print("variants requested:", len(ids))
PY

"$PLINK" --bfile "$BFILE" --extract "$ROOT/data/ld/all759.snplist" \
  --make-just-bim --out "$ROOT/data/ld/all759"
"$PLINK" --bfile "$BFILE" --extract "$ROOT/data/ld/all759.snplist" \
  --r square --out "$ROOT/data/ld/all759ld"

rm -f "$ROOT/data/ld/"*.nosex "$ROOT/data/ld/"*.log "$ROOT/data/ld/all759.snplist"
echo "variants in matrix: $(wc -l < "$ROOT/data/ld/all759.bim")"
echo "done: data/ld/all759.bim, data/ld/all759ld.ld"
