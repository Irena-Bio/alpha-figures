# alpha-figures

Python script for generating alpha-diversity boxplots with Welch’s t-test brackets
for bacterial and eukaryotic microbiome data. It also builds 2×2 composites and
an 8-panel (A–H) figure ready for publication.

## Inputs
Place these Excel files in the repository root (or adjust paths in the script):
- `Bacterial_ASV_alpha.xlsx`
- `Eukaryotic_ASV_alpha.xlsx`

First column = sample IDs (e.g., `KNIN_1`, `KNIN.1.SED`, `KUPC_1`, `VRAB_1`).
Columns containing “ASV” (observed richness) and “shannon” are detected automatically.

## Install
```bash
python -m venv .venv
# Windows: .venv\Scripts\activate
# macOS/Linux:
source .venv/bin/activate
pip install -r requirements.txt

