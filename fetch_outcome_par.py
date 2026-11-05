import csv, os, subprocess, sys, time
from concurrent.futures import ThreadPoolExecutor, as_completed

ROOT = os.environ.get("NPM_MR_ROOT", ".")
N14 = "https://storage.googleapis.com/finngen-public-data-r11/summary_stats/finngen_R11_N14_INFLAMMBREAST.gz"
O15 = "https://storage.googleapis.com/finngen-public-data-r11/summary_stats/finngen_R11_O15_BIRTH_BREAT_INFECT.gz"
DEFAULTS = {"N14": N14, "O15": O15}
if len(sys.argv) < 2:
    sys.exit("usage: fetch_outcome_par.py <N14|O15> [url] [workers]")
name = sys.argv[1]
url = sys.argv[2] if len(sys.argv) > 2 and sys.argv[2].startswith("http") else DEFAULTS[name]
WORKERS = int(sys.argv[3]) if len(sys.argv) > 3 else 12
WINDOW = 250000

pos = []
with open(f"{ROOT}/data/positions_grch38.csv") as fh:
    for r in csv.DictReader(fh):
        pos.append((r["chr"], int(r["pos"]), r["rsid"]))
by_chr = {}
for c, p, rs in pos:
    by_chr.setdefault(c, []).append((p, rs))
for c in by_chr:
    by_chr[c].sort()

regions = []
for c, arr in by_chr.items():
    start = arr[0][0]; prev = arr[0][0]
    for p, _ in arr[1:]:
        if p - prev > WINDOW:
            regions.append((c, start, prev)); start = p
        prev = p
    regions.append((c, start, prev))

target = {rs for _, _, rs in pos}
outf = f"{ROOT}/data/outcome_{name}.csv"

def query(reg):
    c, s, e = reg
    try:
        res = subprocess.run(["tabix", "-s", "1", "-b", "2", "-e", "2", url, f"{c}:{s}-{e}"],
                             capture_output=True, text=True, timeout=120)
        rows = []
        for line in res.stdout.splitlines():
            f = line.split("\t")
            if len(f) < 13: continue
            hit = set(f[4].split(",")) & target
            if hit:
                rows.append(f + [",".join(sorted(hit))])
        return rows
    except Exception:
        return []

allrows = []
t0 = time.time()
with ThreadPoolExecutor(max_workers=WORKERS) as ex:
    futs = {ex.submit(query, r): i for i, r in enumerate(regions)}
    n = 0
    for fu in as_completed(futs):
        allrows.extend(fu.result()); n += 1
        if n % 100 == 0:
            print(n, "/", len(regions), "rows", len(allrows), round(time.time()-t0,1), flush=True)
        if time.time() - t0 > 520 and n > 100:
            break
print("done", n, "regions; rows", len(allrows), flush=True)

with open(outf, "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["chr","pos","ref","alt","rsids","nearest_genes","pval","mlogp","beta","sebeta","af_alt","af_alt_cases","af_alt_controls","query_rsid"])
    w.writerows(allrows)
print("written", len(allrows), flush=True)
