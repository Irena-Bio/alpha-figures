"""
make_alpha_figures.py
Creates group-mean tables and alpha-diversity boxplots with Welch's t-test
brackets for bacteria and eukaryotes, including A–H composite.

Inputs (in the working directory):
  - Bacterial_ASV_alpha.xlsx
  - Eukaryotic_ASV_alpha.xlsx

Outputs (to working directory):
  - *_Alpha_means_by_group.(csv|xlsx)
  - panel_* (bacteria A–D)
  - Euk_* (eukaryotes E–H)
  - Alpha_boxplots_ABCD_with_significance.(png|pdf)
  - Eukaryotic_Alpha_ABCD_with_significance.(png|pdf)
  - Alpha_8panels_A_to_H.(png|pdf)
"""

import re
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from PIL import Image

# ---------- Configuration ----------
IN_BACT = Path("Bacterial_ASV_alpha.xlsx")
IN_EUK  = Path("Eukaryotic_ASV_alpha.xlsx")
OUTDIR  = Path(".")
COLORS  = {"Krka":"#4C78A8", "Kupcina":"#59A14F", "Fish farm":"#F28E2B"}

# Welch t-test (SciPy if available; permutation fallback otherwise)
try:
    from scipy.stats import ttest_ind
    def welch_p(a, b):
        return ttest_ind(a, b, equal_var=False, nan_policy="omit").pvalue
except Exception:
    def welch_p(a, b, n_perms=20000):
        rng = np.random.default_rng(42)
        a = np.asarray(a); b = np.asarray(b)
        a = a[~np.isnan(a)]; b = b[~np.isnan(b)]
        def tstat(x,y):
            return (x.mean()-y.mean())/np.sqrt(x.var(ddof=1)/len(x)+y.var(ddof=1)/len(y))
        t_obs = abs(tstat(a,b))
        pooled = np.concatenate([a,b])
        n = len(a)
        more = 0
        for _ in range(n_perms):
            rng.shuffle(pooled)
            x, y = pooled[:n], pooled[n:]
            if abs(tstat(x,y)) >= t_obs:
                more += 1
        return (more+1)/(n_perms+1)

def stars(p):
    return "***" if p < 0.001 else ("**" if p < 0.01 else ("*" if p < 0.05 else "ns"))

def assign_group(sample: str) -> str:
    """Map sample name to a {River habitat} group.
       Recognizes KNIN/KRKA, KUPC/KUPČ/KUPCINA/KUPČINA, VRAB;
       sediment if endswith 'S' or contains '.SED' (case-insensitive)."""
    su = re.sub(r"\s+", "", str(sample)).upper()
    if any(tok in su for tok in ["KNIN", "KRKA"]):
        river = "Krka"
    elif any(tok in su for tok in ["KUPC", "KUPČ", "KUPCINA", "KUPČINA"]):
        river = "Kupcina"
    elif "VRAB" in su:
        river = "Fish farm"
    else:
        return "Unassigned"
    habitat = "sediment" if (su.endswith("S") or ".SED" in su) else "water"
    return f"{river} {habitat}"

def load_alpha_table(path: Path):
    df = pd.read_excel(path)
    df.columns = [str(c).strip() for c in df.columns]
    sample_col = df.columns[0]
    df[sample_col] = df[sample_col].astype(str).str.strip()
    df["Group"] = df[sample_col].apply(assign_group)

    # Detect metric columns (Observed ASVs & Shannon)
    # Typical names: ASVs, shannon_entropy / shannon
    rich = "ASVs" if "ASVs" in df.columns else next((c for c in df.columns if "asv" in c.lower()), None)
    shan = "shannon_entropy" if "shannon_entropy" in df.columns else (
           "shannon" if "shannon" in df.columns else next((c for c in df.columns if "shannon" in c.lower()), None))
    return df, sample_col, rich, shan

def save_group_means(df: pd.DataFrame, metrics, out_prefix: str):
    wanted = ["Krka water","Krka sediment","Kupcina water","Kupcina sediment","Fish farm water"]
    means = df[df["Group"].isin(wanted)].groupby("Group")[metrics].mean().reindex(wanted)
    means.to_csv(OUTDIR / f"{out_prefix}_Alpha_means_by_group.csv")
    means.to_excel(OUTDIR / f"{out_prefix}_Alpha_means_by_group.xlsx")
    return means

def add_bracket(ax, i, j, y, h, text):
    ax.plot([i, i, j, j], [y, y+h, y+h, y], lw=1.2, c="black")
    ax.text((i+j)/2, y+h*1.05, text, ha="center", va="bottom")

def panel(df, groups, labels, metric, title, ylabel, fname):
    data = [df[df["Group"]==g][metric].dropna().values for g in groups]
    fig, ax = plt.subplots()
    bp = ax.boxplot(data, labels=labels, patch_artist=True)
    for patch, lab in zip(bp['boxes'], labels):
        river = "Fish farm" if lab == "Fish farm" else lab
        patch.set_facecolor(COLORS[river])
        patch.set_alpha(0.9); patch.set_linewidth(1.25)
    for med in bp['medians']: med.set_linewidth(2)

    ax.set_title(title)
    ax.set_ylabel(ylabel)

    # Pairwise Welch brackets
    if any(len(d)==0 for d in data):
        ymax, ymin = 1.0, 0.0
    else:
        ymax = max(np.nanmax(d) for d in data)
        ymin = min(np.nanmin(d) for d in data)
    yr = max(ymax - ymin, 1.0)
    base = ymax + 0.05*yr
    step = 0.06*yr

    pairs = [(i,j) for i in range(1, len(groups)+1) for j in range(i+1, len(groups)+1)]
    pairs = sorted(pairs, key=lambda x: (x[1]-x[0], x[0]))
    for k, (i, j) in enumerate(pairs):
        p = welch_p(data[i-1], data[j-1])
        add_bracket(ax, i, j, base + k*step, step*0.5, stars(p))

    ax.set_ylim(top=base + (len(pairs)+1)*step)
    plt.tight_layout()
    plt.savefig(fname, dpi=300)
    plt.close(fig)

def make_all_panels(alpha_path: Path, out_prefix: str, titles_prefix: str,
                    panel_names=("A","B","C","D")):
    df, sample_col, rich, shan = load_alpha_table(alpha_path)
    metrics = [m for m in [rich, shan] if m is not None]
    save_group_means(df, metrics, out_prefix)

    # Group sets
    sed_groups  = ["Krka sediment","Kupcina sediment"]
    sed_labels  = ["Krka","Kupcina"]
    water_groups = ["Krka water","Kupcina water","Fish farm water"]
    water_labels = ["Krka","Kupcina","Fish farm"]

    # Panel files
    pA = OUTDIR / f"{out_prefix}_{panel_names[0]}_sediment_ASVs_sig.png"
    pB = OUTDIR / f"{out_prefix}_{panel_names[1]}_sediment_Shannon_sig.png"
    pC = OUTDIR / f"{out_prefix}_{panel_names[2]}_water_ASVs_sig.png"
    pD = OUTDIR / f"{out_prefix}_{panel_names[3]}_water_Shannon_sig.png"

    panel(df, sed_groups, sed_labels, rich, f"{panel_names[0]}. {titles_prefix} (Sediment): ASV Richness by Site",
          "Observed ASVs", str(pA))
    panel(df, sed_groups, sed_labels, shan, f"{panel_names[1]}. {titles_prefix} (Sediment): Shannon Diversity by Site",
          "Shannon Index", str(pB))
    panel(df, water_groups, water_labels, rich, f"{panel_names[2]}. {titles_prefix} (Water): ASV Richness by Site",
          "Observed ASVs", str(pC))
    panel(df, water_groups, water_labels, shan, f"{panel_names[3]}. {titles_prefix} (Water): Shannon Diversity by Site",
          "Shannon Index", str(pD))
    return pA, pB, pC, pD

def make_composite_2x2(panels, outfile_png, outfile_pdf):
    imgs = [Image.open(str(p)) for p in panels]
    w = max(i.width for i in imgs); h = max(i.height for i in imgs)
    canvas = Image.new("RGB", (w*2, h*2), "white")
    canvas.paste(imgs[0].resize((w,h)), (0,0))
    canvas.paste(imgs[1].resize((w,h)), (w,0))
    canvas.paste(imgs[2].resize((w,h)), (0,h))
    canvas.paste(imgs[3].resize((w,h)), (w,h))
    canvas.save(outfile_png, format="PNG")
    canvas.save(outfile_pdf, format="PDF")

def make_composite_2x4(bac_panels, euk_panels, outfile_png, outfile_pdf):
    imgs = [Image.open(str(p)) for p in bac_panels + euk_panels]
    w = max(i.width for i in imgs); h = max(i.height for i in imgs)
    canvas = Image.new("RGB", (w*4, h*2), "white")
    for idx, img in enumerate(imgs):
        r = 0 if idx < 4 else 1
        c = idx if idx < 4 else idx - 4
        canvas.paste(img.resize((w,h)), (c*w, r*h))
    canvas.save(outfile_png, format="PNG")
    canvas.save(outfile_pdf, format="PDF")

# ---------- Run ----------
# Bacteria: panels A–D
bac_panels = make_all_panels(IN_BACT, out_prefix="Bacterial", titles_prefix="Bacteria", panel_names=("A","B","C","D"))
make_composite_2x2(bac_panels,
                   OUTDIR / "Alpha_boxplots_ABCD_with_significance.png",
                   OUTDIR / "Alpha_boxplots_ABCD_with_significance.pdf")

# Eukaryotes: panels E–H
euk_panels = make_all_panels(IN_EUK, out_prefix="Eukaryotic", titles_prefix="Eukaryotes", panel_names=("E","F","G","H"))
make_composite_2x2(euk_panels,
                   OUTDIR / "Eukaryotic_Alpha_ABCD_with_significance.png",
                   OUTDIR / "Eukaryotic_Alpha_ABCD_with_significance.pdf")

# 8-panel composite A–H
make_composite_2x4(bac_panels, euk_panels,
                   OUTDIR / "Alpha_8panels_A_to_H.png",
                   OUTDIR / "Alpha_8panels_A_to_H.pdf")
print("Done.")
Add make_alpha_figures.py
