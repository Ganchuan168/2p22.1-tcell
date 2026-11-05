# 2p22.1-tcell

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23954027.svg)](https://doi.org/10.5281/zenodo.23954027)

Analysis code for:

> **A 2p22.1 locus is associated with coordinated T-cell differentiation states and with
> non-puerperal mastitis: a two-sample Mendelian randomisation study of 731 immunophenotypes**

This repository contains the complete analysis pipeline: instrument selection, harmonisation,
forward and reverse two-sample Mendelian randomisation (MR), sensitivity analyses, the
variant-level screen, and generation of every table and figure in the manuscript.

The repository contains **code only**. All input data are public summary statistics
obtained from the GWAS Catalog and FinnGen (Section 3); they are not redistributed here and
are rebuilt by the scripts in Section 5.

---

## 1. Repository layout

```
.
├── scripts/                      analysis code (see Section 5 for run order)
├── data/
│   ├── README.md                 data dictionary: file formats, columns, download commands
│   ├── ivs/                      instrument-level summary statistics (generated; 731 files)
│   ├── reverse/                  reverse-direction associations (generated; 731 files)
│   ├── outcome_*.csv             FinnGen R11 outcome extracts (generated)
│   ├── positions_grch38.csv      GRCh38 coordinates of candidate variants (generated)
│   ├── reverse_instruments.csv   reverse-MR instruments (generated)
│   └── ld/                       output directory for the 1000 Genomes EUR LD matrix
├── results/                      CSV/XLSX analysis output (generated)
├── figures/                      Figure 1–4 source files, PDF and PNG (generated)
└── tables/                       Supplementary Tables S1–S3 (generated)
```

No data files are version-controlled. Every file under `data/` is created by the scripts in
Section 5 from public sources, or is described in `data/README.md`.

---

## 2. Software requirements

| Software | Version used |
|---|---|
| R | 4.6.1 |
| TwoSampleMR | 0.7.11 |
| MRPRESSO | 1.0 |
| ieugwasr | 1.2.0 |
| data.table | 1.18.6.1 |
| dplyr | 1.2.1 |
| ggplot2 | 4.0.3 (with ggrepel 0.9.8, scales 1.4.0, ragg 1.5.2) |
| openxlsx | 4.2.9 |
| Matrix | 1.7.6 |
| PLINK | 1.90b7.2 |
| Python | 3.12 (pandas 2.3.3, numpy 2.5.3, scipy 1.18.1, openpyxl 3.1.5) |
| aria2c, pigz, tabix | for the bulk-download fallback (optional) |

Install the R packages with:

```r
install.packages(c("data.table", "dplyr", "ggplot2", "ggrepel", "scales",
                   "openxlsx", "Matrix", "ragg"))
install.packages("remotes")
remotes::install_github("MRCIEU/TwoSampleMR")
remotes::install_github("rondolab/MR-PRESSO")
```

Python packages are pinned in `requirements.txt` (`pip install -r requirements.txt`).

Package versions are locked in `renv.lock`, which records all 92 packages of the dependency
closure (direct packages and everything they load) at the versions used for the manuscript.
The five GitHub-hosted packages are pinned to the commit they were installed from
(`TwoSampleMR` v0.7.11, `ieugwasr` v1.2.0 and, as hard imports of `TwoSampleMR`, `MR-PRESSO` 1.0,
`MRMix` 0.1.0 and `RadialMR` 1.2.4); no radial-MR or MR-Mix analysis was run, and the scripts
call neither package directly. Each entry carries its commit in `RemoteSha`, so
`renv::restore()` reproduces the analytic stack rather than pulling the current `HEAD`. `sessionInfo.txt` is the complete session recorded with
`Rscript scripts/session_info.R` (R 4.6.1, Ubuntu 24.04, every attached package and every
loaded namespace); re-run that script to refresh it.
`scripts/run_forward.R` uses `parallel::mclapply` on Linux and macOS and falls back to
serial execution on Windows.

---

## 3. Data sources

| Data | Source | Accession / file |
|---|---|---|
| Exposure: 731 circulating immune-cell phenotypes | Orrù et al., *Nat Genet* 2020; 3,757 Sardinian individuals | GWAS Catalog GCST90001391–GCST90002121 |
| Primary outcome: non-puerperal mastitis (ICD-10 N61) | FinnGen release 11 | `N14_INFLAMMBREAST` (2,167 cases / 233,585 controls) |
| Comparator outcome: lactational mastitis (ICD-10 O91) | FinnGen release 11 | `O15_BIRTH_BREAT_INFECT` (1,233 cases / 232,519 controls) |
| Reverse-direction instruments | FinnGen release 5, via OpenGWAS | `finn-b-N14_INFLAMMBREAST` (757 cases / 115,030 controls) |
| LD reference | 1000 Genomes phase 3, European | binary PLINK panel (user-supplied, see below) |

No data files are distributed with this repository. The pipeline downloads the exposure and
outcome summary statistics directly from their public providers (Steps 0 and 2 in Section 5)
and writes them to `data/`, whose layout and column definitions are documented in
`data/README.md`. Downloads from the GWAS Catalog are open; instrument retrieval through
`ieugwasr` needs a free OpenGWAS token (`OPENGWAS_JWT`), and a token-free bulk route over the
GWAS Catalog FTP is provided as an alternative.

---

## 4. Environment variables

All scripts resolve paths relative to a single root directory, so the repository can be
cloned anywhere.

| Variable | Default | Meaning |
|---|---|---|
| `NPM_MR_ROOT` | `.` (current directory) | repository root; **run every command from the repository root** |
| `PLINK_1KG_EUR` | `data/ld/EUR` | prefix of the 1000 Genomes EUR binary panel (`<prefix>.bed/.bim/.fam`) |
| `PLINK_BIN` | `plink1.9` | PLINK 1.9 executable |
| `OPENGWAS_JWT` | — | OpenGWAS token, required only if you re-query instruments with `ieugwasr` |

Example:

```bash
export NPM_MR_ROOT="$PWD"
export PLINK_1KG_EUR=/path/to/1kg_eur/EUR
```

---

## 5. Running the pipeline

Run every command from the repository root, in the order given. Steps 0–2 download the
input data; they must be completed before Step 3.

### Step 0 — instruments (writes `data/ivs/`)

```bash
# Requires OPENGWAS_JWT. Runs in blocks; re-run until all 731 files exist.
Rscript scripts/get_instruments.R 90001391 90002121 3600
```

Alternative bulk route over the GWAS Catalog FTP (no token; needs aria2c and pigz):

```bash
seq 90001391 90002121 | sed 's/^/GCST/' > /tmp/ids.txt
bash scripts/fetch_batch.sh /tmp/ids.txt 8
```

### Step 1 — GRCh38 positions (writes `data/positions_grch38.csv`)

```bash
python3 scripts/get_positions.py
```

### Step 2 — outcome extracts (writes `data/outcome_*.csv`)

```bash
# Re-run until "finished regions" reaches the total number of regions.
python3 scripts/fetch_outcome.py N14 \
  https://storage.googleapis.com/finngen-public-data-r11/summary_stats/finngen_R11_N14_INFLAMMBREAST.gz
python3 scripts/fetch_outcome.py O15 \
  https://storage.googleapis.com/finngen-public-data-r11/summary_stats/finngen_R11_O15_BIRTH_BREAT_INFECT.gz

# or, multi-threaded (writes the same output):
python3 scripts/fetch_outcome_par.py N14 <same url> 12
```

### Step 3 — forward MR

```bash
Rscript scripts/run_forward.R outcome_N14.csv N14
Rscript scripts/run_forward.R outcome_O15.csv O15
```

Instrument selection, LD clumping (r² < 0.001 within 10,000 kb, 1000 Genomes EUR),
MAF filtering (0.01–0.99) and harmonisation all happen inside `build_trait()` in
`scripts/mr_core.R`. With 1 variant a Wald ratio is estimated; with ≥ 2 variants,
inverse-variance weighted (multiplicative random effects) MR plus MR-Egger, weighted
median, weighted mode and simple mode. Benjamini–Hochberg FDR is applied across all
phenotypes with an estimate (`n = 578`) and the Bonferroni threshold is `0.05/578`.

Outputs: `results/forward_N14.csv`, `results/forward_O15.csv`.

### Step 4 — sensitivity analyses

```bash
Rscript scripts/sensitivity.R
```

Cochran's Q, MR-Egger intercept, MR-PRESSO (1,000 simulations), leave-one-out and the
Steiger directionality test for every phenotype with P < 0.05.
Output: `results/sensitivity_hits.csv`.

### Step 5 — reverse-direction MR

```bash
# Requires OPENGWAS_JWT; re-run in blocks. Writes data/reverse/ (731 files).
# The instruments come from the FinnGen release 5 N61 dataset in OpenGWAS
# (finn-b-N14_INFLAMMBREAST), not from release 11; see Section 7.
Rscript scripts/reverse_mr.R 90001391 90002121 3600

Rscript scripts/reverse_analysis.R
```

Output: `results/reverse_NPM.csv`.

### Step 6 — manuscript tables

```bash
Rscript scripts/make_tables.R
```

Outputs: `results/table2.csv`, `results/tables.xlsx` (Tables 1–3).

### Step 7 — figures

```bash
Rscript scripts/fig1_design.R       # Figure 1  study design
Rscript scripts/figures_revised.R   # Figure 2  volcano, Figure 3  forest
Rscript scripts/fig4_locus.R        # Figure 4  2p22.1 locus
```

Outputs: `figures/fig1_design.{pdf,png}`, `figures/fig2_volcano.{pdf,png}`,
`figures/fig3_forest.{pdf,png}`, `figures/fig4_locus.{pdf,png}`, `results/per_variant_wald.csv`.

### Step 8 — variant-level analysis and 2p22.1 LD matrix

```bash
Rscript scripts/make_supp_table.R                 # tables/Supplementary_Table_S1_instrument_level.csv
bash    scripts/make_ld_matrix.sh                 # data/ld/all759.bim, data/ld/all759ld.ld
Rscript scripts/variant_table.R                   # results/variant_level_table.csv
Rscript scripts/correlated_ivw.R                  # results/correlated_ivw_sensitivity.csv
python3 scripts/make_variant_s3.py                # tables/Supplementary_Table_S3_variant_level.csv
```

`make_ld_matrix.sh` needs PLINK 1.9 and a 1000 Genomes phase 3 European binary panel
(`PLINK_1KG_EUR`). Of the 759 variants, 753 are present in that panel; the remaining six
are excluded from the LD-based analyses only.

`variant_table.R` reports the six core variants at 2p22.1 (rs6707421, rs12712610,
rs6751481, rs413431, rs408686, rs17607399), together with their LD cluster (A1, A2, B, C)
and their r with the lead variant rs6707421. `correlated_ivw.R` repeats the IVW estimate
with a generalised-least-squares weighting matrix that accounts for this correlation.
`make_variant_s3.py` builds the variant-level screen: each of the 759 variants is tested
once against N61 through its variant–outcome association, and Benjamini–Hochberg FDR is
applied across the 759 variants (Bonferroni threshold `0.05/759 = 6.59 × 10⁻⁵`).

### Step 9 — supplementary workbook

```bash
python3 scripts/make_supp_workbook.py
```

Assembles Supplementary Tables S1–S3 into a single formatted XLSX workbook.

---

## 6. Output-to-manuscript map

| Manuscript item | Produced by |
|---|---|
| Table 1 (data sources and study parameters) | `scripts/make_tables.R` |
| Table 2 (variant-level association of the six core variants) | `scripts/variant_table.R` |
| Table 3 (phenotype-level ranking, FDR < 0.2) | `scripts/make_tables.R` (`results/table2.csv`) |
| Table 4 (reverse-direction analysis) | `scripts/make_tables.R` |
| Figure 1 (study design) | `scripts/fig1_design.R` |
| Figure 2 (volcano) | `scripts/figures_revised.R` |
| Figure 3 (forest plot) | `scripts/figures_revised.R` |
| Figure 4 (2p22.1 locus) | `scripts/fig4_locus.R` |
| Supplementary Table S1 (instrument level, 1,722 rows) | `scripts/make_supp_table.R` |
| Supplementary Table S2 (phenotype level, 578 rows) | `scripts/make_supp_workbook.py` |
| Supplementary Table S3 (variant level, 759 rows) | `scripts/make_variant_s3.py` |

---

## 7. Notes on the analysis

- **Instrument selection.** Variants were associated with each phenotype at P < 5 × 10⁻⁸,
  clumped at r² < 0.001 within 10,000 kb in 1000 Genomes EUR, and restricted to
  minor allele frequency 0.01–0.99. 612 phenotypes had at least one instrument; 578
  could be harmonised with the outcome and entered the MR analysis.
- **Denominators.** The phenotype-level Bonferroni threshold uses 578 phenotypes
  (`0.05/578 = 8.65 × 10⁻⁵`). Across variants, the 1,722 variant–phenotype pairs
  collapse to 759 distinct variants, so the variant-level threshold is `0.05/759`.
- **Effect scale.** The exposure GWAS inverse-normal transformed the immune traits, so the
  per-allele beta is in standard deviation units and the MR estimate is expressed per one
  standard deviation of the genetically predicted phenotype. The exception is a family of six
  CD45RA/CD28-defined CD8+ T-cell phenotypes (GCST90001695–GCST90001700) analysed on their raw
  measurement scale; five of them entered the analysis, their betas are in the original units,
  and R² is reported as missing for them (`run_mr()` sets R² to NA when the approximation
  exceeds 1). None of the five reached nominal significance. R² for instrument strength is
  calculated on this standardised scale.
- **Reverse-direction instruments.** FinnGen release 11 contains no genome-wide significant
  variant for N61, and the reverse screen also needs the association of each candidate
  variant with all 731 phenotypes. Both requirements are met by the earlier FinnGen release
  held in OpenGWAS (`finn-b-N14_INFLAMMBREAST`, release 5; 757 cases / 115,030 controls),
  from which the seven instruments at P < 5 × 10⁻⁶ were selected. The forward analysis uses
  release 11 throughout.
- **Palindromic variants.** Variants with A/T or C/G alleles and an allele frequency
  between 0.42 and 0.58 were discarded during harmonisation. The frequency used is that of
  the outcome (FinnGen) sample, after any allele flip.
- **Raw-scale phenotypes.** Five phenotypes (two absolute counts and three percentages from
  the CD45RA/CD28-defined CD8+ T-cell family) have exposure betas on their raw measurement
  scale; R² is therefore not reported for those phenotypes in Supplementary Table S1.
- **Random seeds.** MR-PRESSO uses 1,000 simulations; results in the manuscript use the
  default seed of the `MRPRESSO` package.

---

## 8. Data availability and privacy

This repository contains **analysis code only**. It contains:

- no individual-level genotype or phenotype data;
- no participant identifiers, dates of birth, or any other directly identifying field;
- no data derived from the authors' own patients;
- no redistributed summary statistics — every input is downloaded from its public
  provider by the scripts in Section 5.

The exposure data are published, de-identified summary statistics from Orrù et al.
(*Nat Genet* 2020) distributed through the GWAS Catalog; the outcome data are summary
statistics from FinnGen release 11. Both are obtained directly from those providers.

---

## 9. Licence

Code released under the MIT licence (see `LICENSE`).

## 10. Citation

If you use this code, please cite the manuscript and this archive
(`CITATION.cff`), together with the underlying data sources:

- This archive: <https://doi.org/10.5281/zenodo.23954027> (concept DOI: it always resolves to
  the most recent release, so it stays valid across versions)
- Source code: <https://github.com/Ganchuan168/2p22.1-tcell>
- Orrù V, Steri M, Sidore C, et al. Complex genetic signatures in immune cells underlie
  autoimmunity and inform therapy. *Nat Genet* (2020) 52:1036–1045.
- Kurki MI, Karjalainen J, Palta P, et al. FinnGen provides genetic insights from a
  well-phenotyped isolated population. *Nature* (2023) 613:508–518.
- Hemani G, Zheng J, Elsworth B, et al. The MR-Base platform supports systematic causal
  inference across the human phenome. *eLife* (2018) 7:e34408.
