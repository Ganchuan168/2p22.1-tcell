import glob, json, os, time, urllib.request, csv

ROOT = os.environ.get("NPM_MR_ROOT", ".")
CACHE = f"{ROOT}/data/ensembl_cache"
os.makedirs(CACHE, exist_ok=True)

rsids = set()
for f in glob.glob(f"{ROOT}/data/ivs/*.csv"):
    with open(f) as fh:
        for row in csv.DictReader(fh):
            if row.get("rsid"):
                rsids.add(row["rsid"])
rsids = sorted(rsids)
print("unique rsids:", len(rsids), flush=True)

B = 100
for i in range(0, len(rsids), B):
    cf = f"{CACHE}/batch_{i:05d}.json"
    if os.path.exists(cf):
        continue
    chunk = rsids[i:i+B]
    data = json.dumps({"ids": chunk}).encode()
    req = urllib.request.Request("https://rest.ensembl.org/variation/human", data=data,
                                 headers={"Content-Type": "application/json", "Accept": "application/json"})
    d = None
    for attempt in range(4):
        try:
            with urllib.request.urlopen(req, timeout=30) as r:
                d = json.loads(r.read())
            break
        except Exception as e:
            print("retry", i, repr(e)[:80], flush=True)
            time.sleep(2)
    if d is None:
        print("SKIP batch", i, flush=True)
        continue
    with open(cf, "w") as fh:
        json.dump(d, fh)
    print("batch", i, "ok", len(d), flush=True)
    time.sleep(0.3)

out = {}
for cf in glob.glob(f"{CACHE}/batch_*.json"):
    with open(cf) as fh:
        d = json.load(fh)
    for k, v in d.items():
        m = [x for x in v.get("mappings", []) if x.get("assembly_name") == "GRCh38"]
        if m:
            out[k] = {"chr": m[0]["seq_region_name"], "pos": m[0]["start"], "alleles": m[0].get("allele_string", "")}

with open(f"{ROOT}/data/positions_grch38.csv", "w") as fh:
    fh.write("rsid,chr,pos,alleles\n")
    for k, v in sorted(out.items()):
        fh.write(f"{k},{v['chr']},{v['pos']},{v['alleles']}\n")
print("mapped:", len(out), "of", len(rsids), flush=True)
