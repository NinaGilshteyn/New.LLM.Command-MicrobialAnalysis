#!/usr/bin/env python3
"""
run_analysis.py — bundled script for microbiome-analysis skill
==============================================================
Runs complete microbiome analysis: taxonomy → alpha → beta diversity.

Usage:
  python scripts/run_analysis.py \
    --otu otu_table.csv \
    --taxonomy taxonomy.csv \
    --metadata metadata.csv \
    --group-col group \
    --output-dir results/
"""

import argparse
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scipy.spatial.distance import braycurtis
from scipy.stats import mannwhitneyu, entropy
from sklearn.manifold import MDS
from pathlib import Path
import warnings
warnings.filterwarnings("ignore")

PHYLUM_COLORS = {
    "Firmicutes":      "#4A90D9",
    "Bacteroidetes":   "#50C878",
    "Actinobacteria":  "#F5A623",
    "Proteobacteria":  "#E84855",
    "Verrucomicrobia": "#9B59B6",
    "Fusobacteria":    "#95A5A6",
    "Tenericutes":     "#BDC3C7",
    "Spirochaetes":    "#D5DBDB",
}
GROUP_PALETTE = ["#2E86AB", "#E84855", "#3BB273", "#F18F01", "#7B2D8B"]


def load_data(otu_path, tax_path, meta_path, group_col):
    otu  = pd.read_csv(otu_path,  index_col=0)
    tax  = pd.read_csv(tax_path,  index_col=0)
    meta = pd.read_csv(meta_path, index_col=0)
    # Align: keep only samples present in all three files
    shared = otu.index.intersection(meta.index)
    otu  = otu.loc[shared]
    meta = meta.loc[shared]
    groups = meta[group_col].unique().tolist()
    return otu, tax, meta, group_col, groups


def relative_abundance(otu):
    return otu.div(otu.sum(axis=1), axis=0)


def collapse_to_phylum(rel, tax):
    df = rel.T.copy()
    df["phylum"] = tax["phylum"]
    return df.groupby("phylum").sum().T


def shannon(row):
    p = row[row > 0]
    return float(entropy(p))

def richness(row):
    return int((row > 0).sum())

def compute_alpha(rel, meta, group_col):
    alpha = pd.DataFrame({
        "shannon":  rel.apply(shannon, axis=1),
        "richness": rel.apply(richness, axis=1),
        "group":    meta[group_col],
    })
    return alpha

def compute_beta(rel):
    samples = rel.index.tolist()
    n = len(samples)
    dm = np.zeros((n, n))
    for i in range(n):
        for j in range(i+1, n):
            d = braycurtis(rel.iloc[i].values, rel.iloc[j].values)
            dm[i,j] = d
            dm[j,i] = d
    mds = MDS(n_components=2, dissimilarity="precomputed", random_state=42, normalized_stress=False)
    coords = mds.fit_transform(dm)
    coords_df = pd.DataFrame(coords, index=samples, columns=["PC1","PC2"])
    return dm, coords_df


def make_figure(phylum_abund, alpha, coords_df, meta, group_col, groups, outpath):
    color_map = {g: GROUP_PALETTE[i % len(GROUP_PALETTE)] for i, g in enumerate(groups)}
    fig, axes = plt.subplots(1, 3, figsize=(15, 5.5))
    fig.patch.set_facecolor("white")

    # ── Panel 1: taxonomy stacked bar ──
    ax1 = axes[0]
    sample_order = []
    boundaries = []
    for g in groups:
        g_samples = meta[meta[group_col] == g].index.tolist()
        boundaries.append((len(sample_order), len(sample_order) + len(g_samples), g))
        sample_order.extend(g_samples)

    plot_data = phylum_abund.loc[sample_order]
    bottom = np.zeros(len(plot_data))
    phyla_order = plot_data.mean().sort_values(ascending=False).index

    for phylum in phyla_order:
        vals = plot_data[phylum].values
        color = PHYLUM_COLORS.get(phylum, "#CCCCCC")
        ax1.bar(range(len(plot_data)), vals, bottom=bottom,
                color=color, label=phylum, width=0.85, linewidth=0)
        bottom += vals

    for start, end, g in boundaries:
        mid = (start + end - 1) / 2
        ax1.axvline(x=end - 0.5, color="black", linewidth=0.8, linestyle="--", alpha=0.4)
        ax1.text(mid, 1.03, g, ha="center", fontsize=8,
                 color=color_map[g], fontweight="bold")

    ax1.set_xticks(range(len(sample_order)))
    ax1.set_xticklabels(sample_order, rotation=45, ha="right", fontsize=7)
    ax1.set_ylabel("Relative abundance", fontsize=10)
    ax1.set_title("Phylum-level composition", fontsize=11, fontweight="bold")
    ax1.set_ylim(0, 1.12)
    ax1.legend(loc="upper right", fontsize=7, framealpha=0.8,
               title="Phylum", title_fontsize=8)
    ax1.spines[["top","right"]].set_visible(False)

    # ── Panel 2: alpha diversity ──
    ax2 = axes[1]
    positions = list(range(1, len(groups) + 1))
    pvals = []

    for pos, g in zip(positions, groups):
        vals = alpha[alpha["group"] == g]["shannon"].values
        bp = ax2.boxplot(vals, positions=[pos], widths=0.45, patch_artist=True,
                         medianprops=dict(color="white", linewidth=2),
                         whiskerprops=dict(color=color_map[g], linewidth=1.2),
                         capprops=dict(color=color_map[g], linewidth=1.2),
                         flierprops=dict(marker="o", color=color_map[g], alpha=0.5, ms=5))
        bp["boxes"][0].set_facecolor(color_map[g])
        bp["boxes"][0].set_alpha(0.7)
        jitter = np.random.uniform(-0.1, 0.1, len(vals))
        ax2.scatter(np.full(len(vals), pos) + jitter, vals,
                    color=color_map[g], alpha=0.8, s=28, zorder=3)

    if len(groups) == 2:
        stat, pval = mannwhitneyu(
            alpha[alpha["group"] == groups[0]]["shannon"].values,
            alpha[alpha["group"] == groups[1]]["shannon"].values,
            alternative="two-sided"
        )
        y_max = alpha["shannon"].max() + 0.25
        ax2.plot([1,1,2,2],[y_max, y_max+0.05, y_max+0.05, y_max], color="black", lw=1)
        sig = "***" if pval<0.001 else "**" if pval<0.01 else "*" if pval<0.05 else "ns"
        ax2.text(1.5, y_max+0.08, f"{sig}\np={pval:.4f}", ha="center", fontsize=9)

    ax2.set_xticks(positions)
    ax2.set_xticklabels(groups, fontsize=10)
    ax2.set_ylabel("Shannon entropy (H′)", fontsize=10)
    ax2.set_title("Alpha diversity", fontsize=11, fontweight="bold")
    ax2.spines[["top","right"]].set_visible(False)

    # ── Panel 3: beta diversity PCoA ──
    ax3 = axes[2]
    coords_df["group"] = meta[group_col]
    for g in groups:
        mask = coords_df["group"] == g
        ax3.scatter(coords_df.loc[mask,"PC1"], coords_df.loc[mask,"PC2"],
                    c=color_map[g], s=70, alpha=0.85,
                    edgecolors="white", linewidths=0.6, label=g, zorder=3)

    patches = [mpatches.Patch(color=color_map[g], label=g, alpha=0.8) for g in groups]
    ax3.legend(handles=patches, fontsize=9, framealpha=0.8)
    ax3.set_xlabel("PCoA Axis 1", fontsize=10)
    ax3.set_ylabel("PCoA Axis 2", fontsize=10)
    ax3.set_title("Beta diversity — Bray-Curtis PCoA", fontsize=11, fontweight="bold")
    ax3.axhline(0, color="#dddddd", lw=0.5)
    ax3.axvline(0, color="#dddddd", lw=0.5)
    ax3.spines[["top","right"]].set_visible(False)

    plt.suptitle(f"Microbiome analysis: {' vs '.join(groups)}",
                 fontsize=13, fontweight="bold", y=1.01)
    plt.tight_layout()
    plt.savefig(outpath, dpi=180, bbox_inches="tight")
    print(f"✓ Figure saved: {outpath}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--otu",       required=True)
    parser.add_argument("--taxonomy",  required=True)
    parser.add_argument("--metadata",  required=True)
    parser.add_argument("--group-col", default="group")
    parser.add_argument("--output-dir",default="results")
    args = parser.parse_args()

    out = Path(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)

    otu, tax, meta, group_col, groups = load_data(
        args.otu, args.taxonomy, args.metadata, args.group_col
    )
    print(f"Loaded: {otu.shape[0]} samples × {otu.shape[1]} OTUs | Groups: {groups}")

    rel         = relative_abundance(otu)
    phylum      = collapse_to_phylum(rel, tax)
    alpha       = compute_alpha(rel, meta, group_col)
    dm, coords  = compute_beta(rel)

    # Save CSVs
    alpha.to_csv(out / "alpha_diversity.csv")
    pd.DataFrame(dm, index=otu.index, columns=otu.index).to_csv(out / "beta_diversity_matrix.csv")
    phylum.to_csv(out / "phylum_abundance.csv")

    # Save figure
    make_figure(phylum, alpha, coords, meta, group_col, groups,
                out / "microbiome_analysis.png")

    # Print summary
    print("\n=== Alpha diversity summary ===")
    print(alpha.groupby("group")[["shannon","richness"]].mean().round(3))
    print("\n=== Outputs ===")
    for f in sorted(out.iterdir()):
        print(f"  {f}")


if __name__ == "__main__":
    main()
