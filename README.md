# GO_Analysis

Standalone GO analysis + visualization workflow for Drosophila gene sets.

## Included scripts
- `GenerateGOChartReport.py`: runs DAVID GO enrichment on all CSV/XLSX files in an input folder and writes GO output tables.
- `VisualizeGOResults.py`: creates bubble plots from GO output tables (default mode and `--by_GO` mode).
- `ProcessGOresults.py`: post-processes GO result workbooks written by the generator.

## Included runtime data
- `Data/fb_synonym_fb_2024_06.tsv`
- `Data/GO-Term_Category-Mapping.xlsx` (default category/mapping file)
- `Data/GO-Term_Category-Mapping_Corr_Sets.xlsx` (alternate mapping file for correlation sets)
- `Data/GO-TEST/*.csv` (test gene sets)

## Requirements
- Python 3.10+
- Python packages: `pandas`, `numpy`, `matplotlib`, `openpyxl`, `suds`
- DAVID web service access with a registered email

## Quick start
From the repo root:

```bash
python3 GenerateGOChartReport.py Data/GO-TEST 0
python3 VisualizeGOResults.py Data/GO-TEST/GO_Analysis
python3 VisualizeGOResults.py Data/GO-TEST/GO_Analysis --by_GO
```

Optional DAVID account override:

```bash
python3 GenerateGOChartReport.py Data/GO-TEST 0 --david-email "you@example.org"
```

## Outputs
- GO result workbooks: `Data/GO-TEST/GO_Analysis/*.xlsx`
- Plots:
  - Default mode: `Data/GO-TEST/GO_Analysis/Plots/*.png`
  - `--by_GO` mode: `Data/GO-TEST/GO_Analysis/Plots/<file_stem>/All_Aspects.png` and per-aspect PNGs
