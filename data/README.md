# Data dictionary

This directory holds the input data used by the pipeline. Nothing here is
version-controlled: every file is written by the scripts in Section 5 of the main
`README.md`, from public, aggregate summary statistics. No individual-level data are
used at any step, and no file contains a participant identifier of any kind.

## `ivs/` — instrument-level exposure statistics (731 files)

One CSV per exposure phenotype, named by its GWAS Catalog study accession
(`ebi-a-GCST90001391.csv` … `ebi-a-GCST90002121.csv`). Written by
`scripts/get_instruments.R`; the released analysis used one row per instrument–phenotype
pair (2,814 rows, 1,561 distinct variants).

| Column | Meaning |
|---|---|
| `id` | GWAS Catalog study accession |
| `trait` | phenotype name |
| `chr`, `position` | GRCh37 coordinate |
| `rsid` | variant identifier |
| `ea`, `nea` | effect allele, other allele |
| `eaf` | effect-allele frequency |
| `beta`, `se` | effect size and standard error per allele |
| `p` | P value for the variant–phenotype association |
| `n` | sample size (3,658 for the 731 phenotypes) |

Only variants reaching P < 5 × 10⁻⁸ are retained, so each file is small.

**Source.** Orrù V, Steri M, Sidore C, et al. *Nat Genet* (2020) 52:1036–1045;
GWAS Catalog accessions GCST90001391–GCST90002121.
**Regenerate with.** `Rscript scripts/get_instruments.R 90001391 90002121 <seconds>`
(needs `OPENGWAS_JWT`), or the token-free bulk route
`bash scripts/fetch_batch.sh <id_list> 8` followed by the same filtering.

## `reverse/` — reverse-direction associations (731 files)

Associations between the 7 NPM-associated variants (`reverse_instruments.csv`) and each
of the 731 phenotypes, retrieved from OpenGWAS. The N61 source dataset is
`finn-b-N14_INFLAMMBREAST`, which corresponds to **FinnGen release 5** (757 cases /
115,030 controls), not release 11: release 11 was used for the forward analysis, but the
reverse screen also needs the association of each candidate variant with all 731
phenotypes, which that resource provides in one place. Columns follow the
`ieugwasr::associations()` output; `proxy = TRUE` marks rows where a correlated proxy was
returned in place of the requested variant. `reverse_analysis.R` discards proxy rows, as
the forward analysis uses no proxies.

**Regenerate with.** `Rscript scripts/reverse_mr.R 90001391 90002121 <seconds>`.

## `outcome_N14.csv` — FinnGen R11 `N14_INFLAMMBREAST` extract

Non-puerperal mastitis (ICD-10 N61): 2,167 cases / 233,585 controls.

## `outcome_O15.csv` — FinnGen R11 `O15_BIRTH_BREAT_INFECT` extract

Lactational mastitis (ICD-10 O91): 1,233 cases / 232,519 controls.

Both outcome files contain every FinnGen variant within ±250 kb of a candidate exposure
variant, so that harmonisation can look up each instrument by position as well as by rsID.

| Column | Meaning |
|---|---|
| `chr`, `pos`, `ref`, `alt` | GRCh38 coordinate and alleles |
| `rsids` | rsID(s) at the position |
| `nearest_genes` | nearest gene(s) |
| `pval`, `mlogp` | P value and −log₁₀(P) |
| `beta`, `sebeta` | effect size and standard error for `alt` |
| `af_alt` | `alt` allele frequency in the combined sample |
| `af_alt_cases`, `af_alt_controls` | `alt` frequency in cases and controls |
| `query_rsid` | the candidate variant this row was retrieved for |

**Source.** FinnGen release 11 (<https://www.finngen.fi/>).
**Regenerate with.** `python3 scripts/fetch_outcome.py N14 <url>` or
`python3 scripts/fetch_outcome_par.py N14 <url> 12`.

## `positions_grch38.csv`

GRCh38 coordinates and allele strings for all 1,561 candidate variants, obtained from the
Ensembl REST API (`/variation/human`). Used to place variants on the genome before the
outcome lookup.
**Regenerate with.** `python3 scripts/get_positions.py`.

## `reverse_instruments.csv`

The 7 variants associated with NPM at P < 5 × 10⁻⁶ after LD clumping, with their
NPM effect size and standard error. Input to `scripts/reverse_analysis.R`.
**Regenerate with.** `Rscript scripts/reverse_mr.R 90001391 90002121 <seconds>`.

## `ld/`

Output directory for the 1000 Genomes phase 3 European LD matrix
(`all759.bim`, `all759ld.ld`), built by `scripts/make_ld_matrix.sh`. The matrix itself is
not version-controlled because it is large and fully reproducible from the public
reference panel.
