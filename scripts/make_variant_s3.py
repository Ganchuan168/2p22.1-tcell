"""Build Supplementary Table S3: one row per variant used in the variant-level
screen (759 variants).  Each variant is tested once against N61: the Wald-ratio
P value equals the variant-outcome P value, identical across every phenotype the
variant instruments, so no phenotype selection enters the test."""
import numpy as np
import pandas as pd
from scipy.stats import norm
import os

ROOT = os.environ.get("NPM_MR_ROOT", ".")
s1 = pd.read_csv(f"{ROOT}/tables/Supplementary_Table_S1_instrument_level.csv")

# exposure row for each variant: strongest exposure association
s1 = s1.sort_values("p_exposure", kind="mergesort")
first = s1.drop_duplicates("rsid", keep="first").copy()
n_phen = s1.groupby("rsid")["trait_id"].nunique().rename("n_phenotypes")

v = first.merge(n_phen, on="rsid", how="left")
v["p_out"] = 2 * norm.sf(np.abs(v["beta_outcome"] / v["se_outcome"]))
v["OR"] = np.exp(v["beta_outcome"])
v["L95"] = np.exp(v["beta_outcome"] - 1.959964 * v["se_outcome"])
v["U95"] = np.exp(v["beta_outcome"] + 1.959964 * v["se_outcome"])

# BH across the 759 variants
order = np.argsort(v["p_out"].values)
m = len(v)
q = np.empty(m)
prev = 1.0
for rank, idx in enumerate(order[::-1], start=1):
    k = m - rank + 1
    val = min(prev, v["p_out"].values[idx] * m / k)
    q[idx] = val
    prev = val
v["fdr"] = q

cluster = {
    "rs6707421": "A1", "rs12712610": "A1",
    "rs6751481": "A2",
    "rs413431": "B", "rs408686": "B",
    "rs17607399": "C",
}
v["cluster"] = v["rsid"].map(cluster).fillna("")

out = v[["rsid", "chr", "position_grch38", "effect_allele", "other_allele",
         "eaf", "beta_exposure", "se_exposure", "beta_outcome", "se_outcome",
         "OR", "L95", "U95", "p_out", "fdr", "cluster", "n_phenotypes", "phenotype"]].copy()
out.columns = ["Variant", "Chromosome", "Position (GRCh38)", "EA", "OA", "EAF",
               "Beta exposure", "SE exposure", "Beta N61", "SE N61",
               "OR N61", "95% CI lower", "95% CI upper", "P N61", "FDR (BH)",
               "LD cluster", "Phenotypes instrumented", "Exposure phenotype"]
# For the six core variants, use the exposure phenotype and effects reported in Table 2
core = pd.read_csv(f"{ROOT}/results/variant_level_table.csv")
core = core.set_index("variant")
for rs, row in core.iterrows():
    m = out["Variant"] == rs
    if m.any():
        out.loc[m, "EA"] = row["ea"]
        out.loc[m, "OA"] = row["oa"]
        out.loc[m, "EAF"] = row["eaf"]
        out.loc[m, "Beta exposure"] = row["b_exp"]
        out.loc[m, "SE exposure"] = row["se_exp"]
        out.loc[m, "Exposure phenotype"] = row["phenotype"]

out = out.sort_values(["Chromosome", "Position (GRCh38)"])
out.to_csv(f"{ROOT}/tables/Supplementary_Table_S3_variant_level.csv", index=False)
print("variants", len(out), "Bonferroni hits", int((out["P N61"] < 0.05 / len(out)).sum()),
      "FDR<0.05", int((out["FDR (BH)"] < 0.05).sum()))
