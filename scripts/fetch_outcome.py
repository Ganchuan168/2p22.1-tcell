import csv, os, subprocess, sys, time

ROOT = os.environ.get("NPM_MR_ROOT", ".")
N14 = "https://storage.googleapis.com/finngen-public-data-r11/summary_stats/finngen_R11_N14_INFLAMMBREAST.gz"
O15 = "https://storage.googleapis.com/finngen-public-data-r11/summary_stats/finngen_R11_O15_BIRTH_BREAT_INFECT.gz"

name, url = sys.argv[1], sys.argv[2]
outf = f"{ROOT}/data/outcome_{name}.csv"
donef = f"{ROOT}/data/outcome_{name}.done"
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
print("regions:", len(regions), flush=True)

target = {rs for _, _, rs in pos}
seen = set()
done_regions = 0
if os.path.exists(donef):
    done_regions = int(open(donef).read().strip())
newfile = not os.path.exists(outf)
fh_out = open(outf, "a", newline="")
w = csv.writer(fh_out)
if newfile:
    w.writerow(["chr", "pos", "ref", "alt", "rsids", "nearest_genes", "pval", "mlogp", "beta", "sebeta", "af_alt", "af_alt_cases", "af_alt_controls", "query_rsid"])

t0 = time.time()
for i, (c, s, e) in enumerate(regions):
    if i < done_regions:
        continue
    if time.time() - t0 > 500:
        print("budget stop at region", i, flush=True)
        break
    try:
        res = subprocess.run(["tabix", "-s", "1", "-b", "2", "-e", "2", url, f"{c}:{s}-{e}"],
                             capture_output=True, text=True, timeout=90)
    except subprocess.TimeoutExpired:
        print("timeout region", i, c, s, e, flush=True)
        continue
    for line in res.stdout.splitlines():
        f = line.split("\t")
        if len(f) < 13:
            continue
        ids = set(f[4].split(","))
        hit = ids & target
        if hit:
            w.writerow(f + [",".join(sorted(hit))])
            seen |= hit
    done_regions = i + 1
    open(donef, "w").write(str(done_regions))
    if i % 50 == 0:
        print("region", i, "/", len(regions), "matched", len(seen), flush=True)

fh_out.close()
print("finished regions", done_regions, "matched rsids", len(seen), "of", len(target), flush=True)
