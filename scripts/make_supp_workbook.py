import pandas as pd, numpy as np
from openpyxl import Workbook
from openpyxl.styles import Font, Alignment, PatternFill, Border, Side
from openpyxl.utils import get_column_letter
from openpyxl.utils.dataframe import dataframe_to_rows
import os

ROOT = os.environ.get("NPM_MR_ROOT", ".")

s1 = pd.read_csv(f"{ROOT}/tables/Supplementary_Table_S1_instrument_level.csv")
s1.loc[s1["R2"] > 1, "R2"] = np.nan
f = pd.read_csv(f"{ROOT}/results/forward_N14.csv")
sen = pd.read_csv(f"{ROOT}/results/sensitivity_hits.csv")
s3 = pd.read_csv(f"{ROOT}/tables/Supplementary_Table_S3_variant_level.csv")

f = f[f["p"].notna()].copy()
s2 = f.merge(sen[["id", "presso_global_p", "presso_outliers", "loo_max_z", "steiger_dir"]],
             on="id", how="left")
s2 = s2.sort_values(["p", "trait"]).reset_index(drop=True)
s2.loc[s2["R2"] > 1, "R2"] = np.nan
s2 = s2[["trait", "id", "nsnp", "Fmean", "Fmin", "R2", "method", "b", "se", "p", "OR", "L95", "U95",
         "FDR", "direction", "egger_intercept", "egger_intercept_p", "Q_p",
         "egger_b", "egger_p", "wm_b", "wm_p", "wmode_b", "wmode_p", "smode_b", "smode_p",
         "presso_global_p", "presso_outliers", "loo_max_z", "steiger_dir"]]
s2.columns = ["Phenotype", "GWAS trait ID", "No. of variants", "Mean F", "Minimum F", "Cumulative R2",
              "Estimator", "Beta", "SE", "P", "OR", "95% CI lower", "95% CI upper", "FDR (BH)",
              "Direction", "MR-Egger intercept", "Egger intercept P", "Cochran Q P",
              "MR-Egger beta", "MR-Egger P", "Weighted median beta", "Weighted median P",
              "Weighted mode beta", "Weighted mode P", "Simple mode beta", "Simple mode P",
              "MR-PRESSO global P", "MR-PRESSO outliers", "LOO max |z|", "Steiger direction"]

s1.columns = ["Phenotype", "GWAS trait ID", "Variant", "Chromosome", "Position (GRCh38)",
              "Effect allele", "Other allele", "EAF", "Beta (exposure)", "SE (exposure)",
              "P (exposure)", "Beta (outcome)", "SE (outcome)", "F statistic", "R2"]

wb = Workbook()
thin = Side(style="thin", color="D9D9D9")
hdr_fill = PatternFill("solid", fgColor="1F3864")
hdr_font = Font(bold=True, color="FFFFFF", size=10)
body_font = Font(size=10)

def add_sheet(name, df, widths, title, fcols=()):
    ws = wb.create_sheet(name)
    ws["A1"] = title
    ws["A1"].font = Font(bold=True, size=12)
    ws.merge_cells(start_row=1, start_column=1, end_row=1, end_column=len(df.columns))
    ws.row_dimensions[1].height = 22
    for j, c in enumerate(df.columns, 1):
        cell = ws.cell(row=3, column=j, value=c)
        cell.font = hdr_font; cell.fill = hdr_fill
        cell.alignment = Alignment(horizontal="center", vertical="center", wrap_text=True)
        ws.column_dimensions[get_column_letter(j)].width = widths[j-1]
    ws.row_dimensions[3].height = 34
    for i, row in enumerate(df.itertuples(index=False), 4):
        for j, v in enumerate(row, 1):
            if isinstance(v, (float, np.floating)) and np.isnan(v):
                v = None
            cell = ws.cell(row=i, column=j, value=v)
            cell.font = body_font
            cell.border = Border(bottom=thin)
            if j in fcols and v is not None:
                cell.number_format = "0.00E+00"
    ws.freeze_panes = "C4"
    ws.auto_filter.ref = f"A3:{get_column_letter(len(df.columns))}{3+len(df)}"
    return ws

# S1
ws1 = add_sheet("Table S1", s1,
    [46, 20, 12, 8, 14, 10, 10, 8, 13, 13, 13, 13, 13, 10, 10],
    "Supplementary Table S1. Instrument-level data for all 1,722 variants used in the forward analysis",
    fcols=(9, 10, 11, 12, 13))
for i in range(4, 4+len(s1)):
    ws1.cell(row=i, column=11).number_format = "0.00E+00"
    ws1.cell(row=i, column=8).number_format = "0.000"

# S2
ws2 = add_sheet("Table S2", s2,
    [46, 20, 8, 9, 9, 13, 12, 11, 10, 11, 9, 11, 11, 11, 10, 12, 12, 11,
     11, 11, 12, 12, 12, 12, 11, 11, 12, 12, 10, 10],
    "Supplementary Table S2. Phenotype-level Mendelian randomization estimates and sensitivity statistics (578 phenotypes)",
    fcols=(10, 14, 17, 18, 20, 22, 24, 26, 27))
for i in range(4, 4+len(s2)):
    for c in (8, 9, 12, 13, 16, 19, 21, 23, 25):
        ws2.cell(row=i, column=c).number_format = "0.000"
    for c in (10, 14, 17, 18, 20, 22, 24, 26, 27):
        ws2.cell(row=i, column=c).number_format = "0.00E+00"

# S3
s3 = s3.sort_values(["P N61"]).reset_index(drop=True)
ws3 = add_sheet("Table S3", s3,
    [12, 10, 13, 6, 6, 8, 11, 11, 11, 10, 10, 10, 10, 11, 10, 9, 11, 46],
    "Supplementary Table S3. Variant-level association with N61-coded non-puerperal mastitis for all 759 variants used in the screen (ordered by P)",
    fcols=(14, 15))
for i in range(4, 4+len(s3)):
    ws3.cell(row=i, column=6).number_format = "0.000"
    for c in (7, 8, 9, 10):
        ws3.cell(row=i, column=c).number_format = "0.0000"
    for c in (11, 12, 13):
        ws3.cell(row=i, column=c).number_format = "0.000"
    for c in (14, 15):
        ws3.cell(row=i, column=c).number_format = "0.00E+00"

# Notes
ws = wb.create_sheet("Notes")
ws.column_dimensions["A"].width = 24
ws.column_dimensions["B"].width = 118
notes = [
 ("Supplementary Table S1", "One row per variant-phenotype pair. Beta and SE for the exposure are on the standardized (per-SD) scale reported by the GWAS Catalog; beta and SE for the outcome are log odds ratios per allele after harmonisation to the exposure effect allele. F = (beta/SE)^2 for the exposure; R2 = 2 x EAF x (1 - EAF) x beta^2."),
 ("Supplementary Table S2", "One row per phenotype analysed (578). OR and 95% CI are exp(beta) from the inverse-variance weighted estimator (Wald ratio for single-variant phenotypes). FDR is the Benjamini-Hochberg adjusted P across the 578 phenotypes. Sensitivity statistics are reported where estimable: MR-Egger and Cochran's Q require at least three variants, MR-PRESSO at least four, and leave-one-out at least two. MR-PRESSO, leave-one-out and Steiger statistics were computed for the 30 phenotypes with P < 0.05."),
 ("Supplementary Table S3", "One row per variant used in the variant-level screen (759 variants). Beta N61 and SE N61 are the per-allele log odds ratio and its standard error from the variant-outcome association after harmonisation to the exposure effect allele; OR = exp(Beta N61). For a single-variant phenotype the Wald-ratio P value equals this variant-outcome P value and is identical across every phenotype the variant instruments, so each variant is tested once and no phenotype selection enters the test. FDR is the Benjamini-Hochberg value across the 759 variants; the Bonferroni threshold is 0.05/759 = 6.59 x 10^-5. LD cluster refers to the four clusters at 2p22.1 defined in the manuscript; blank cells indicate variants outside that region."),
 ("Thresholds", "Phenotype-level Bonferroni threshold 0.05/578 = 8.65 x 10^-5. Variant-level Bonferroni threshold 0.05/759 = 6.59 x 10^-5."),
 ("R2", "R2 is not reported for five absolute-count phenotypes because their exposure betas were on a raw count scale and the approximation 2 x EAF x (1 - EAF) x beta^2 exceeds 1. The F statistic is unaffected because it is scale-free. Summary statistics in the manuscript use the 573 phenotypes with an estimable R2."),
 ("Abbreviations", "AC, absolute count; EAF, effect allele frequency; FDR, false discovery rate; LOO, leave-one-out; MFI, median fluorescence intensity; MR, Mendelian randomization; OR, odds ratio; RC, relative count; SD, standard deviation."),
]
ws["A1"] = "Supplementary tables S1-S3: column definitions and notes"
ws["A1"].font = Font(bold=True, size=12)
ws.merge_cells("A1:B1")
for i, (a, b) in enumerate(notes, 3):
    ws.cell(row=i, column=1, value=a).font = Font(bold=True, size=10)
    ws.cell(row=i, column=1).alignment = Alignment(vertical="top")
    c = ws.cell(row=i, column=2, value=b)
    c.font = Font(size=10); c.alignment = Alignment(wrap_text=True, vertical="top")
    ws.row_dimensions[i].height = 15 * (len(b)//105 + 1) + 10

del wb["Sheet"]
out = f"{ROOT}/tables/Supplementary_Tables_S1_S3.xlsx"
wb.save(out)
print("saved", out, len(s1), len(s2), len(s3))
