#!/usr/bin/env python3
"""
generate_synthetic.py — bundled script for microbiome-analysis skill
=====================================================================
Generates a realistic synthetic OTU table + taxonomy + metadata.

Usage:
  python scripts/generate_synthetic.py \
    --n-samples 20 \
    --n-otus 20 \
    --groups "Healthy,Antibiotic" \
    --output-dir .
"""

import argparse
import numpy as np
import pandas as pd
from pathlib import Path

TAXA_LIBRARY = [
    ("OTU_001","Faecalibacterium","Ruminococcaceae","Firmicutes",    180),
    ("OTU_002","Ruminococcus",    "Ruminococcaceae","Firmicutes",    120),
    ("OTU_003","Blautia",         "Lachnospiraceae","Firmicutes",    100),
    ("OTU_004","Roseburia",       "Lachnospiraceae","Firmicutes",     90),
    ("OTU_005","Lactobacillus",   "Lactobacillaceae","Firmicutes",    30),
    ("OTU_006","Streptococcus",   "Streptococcaceae","Firmicutes",    25),
    ("OTU_007","Clostridium",     "Clostridiaceae","Firmicutes",      40),
    ("OTU_008","Eubacterium",     "Eubacteriaceae","Firmicutes",      60),
    ("OTU_009","Bacteroides",     "Bacteroidaceae","Bacteroidetes",  200),
    ("OTU_010","Prevotella",      "Prevotellaceae","Bacteroidetes",  150),
    ("OTU_011","Parabacteroides", "Tannerellaceae","Bacteroidetes",   60),
    ("OTU_012","Alistipes",       "Rikenellaceae","Bacteroidetes",    50),
    ("OTU_013","Bifidobacterium", "Bifidobacteriaceae","Actinobacteria",40),
    ("OTU_014","Collinsella",     "Eggerthellaceae","Actinobacteria", 30),
    ("OTU_015","Escherichia",     "Enterobacteriaceae","Proteobacteria",5),
    ("OTU_016","Klebsiella",      "Enterobacteriaceae","Proteobacteria",3),
    ("OTU_017","Akkermansia",     "Akkermansiaceae","Verrucomicrobia",35),
    ("OTU_018","Fusobacterium",   "Fusobacteriaceae","Fusobacteria",   2),
    ("OTU_019","Mycoplasma",      "Mycoplasmataceae","Tenericutes",    1),
    ("OTU_020","Treponema",       "Spirochaetaceae","Spirochaetes",    1),
]

# Group-specific abundance modifiers (multipliers on baseline)
GROUP_MODIFIERS = {
    "Healthy":    {"Firmicutes":1.0,"Bacteroidetes":1.0,"Proteobacteria":0.3,"Actinobacteria":1.0},
    "Antibiotic": {"Firmicutes":0.4,"Bacteroidetes":0.3,"Proteobacteria":8.0,"Actinobacteria":0.6},
    "Disease":    {"Firmicutes":0.6,"Bacteroidetes":0.7,"Proteobacteria":3.0,"Actinobacteria":0.8},
    "Control":    {"Firmicutes":1.0,"Bacteroidetes":1.0,"Proteobacteria":0.3,"Actinobacteria":1.0},
}
DEFAULT_MODIFIER = {"Firmicutes":1.0,"Bacteroidetes":1.0,"Proteobacteria":1.0,"Actinobacteria":1.0}


def simulate_group(mean_vec, n_samples, noise=0.35):
    counts = []
    for _ in range(n_samples):
        c = np.random.negative_binomial(
            n=max(1/noise, 0.1),
            p=1/(1 + mean_vec * noise),
        ).astype(float)
        counts.append(c)
    return np.array(counts)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--n-samples", type=int, default=20)
    parser.add_argument("--n-otus",    type=int, default=20)
    parser.add_argument("--groups",    type=str, default="Healthy,Antibiotic")
    parser.add_argument("--output-dir",type=str, default=".")
    parser.add_argument("--seed",      type=int, default=42)
    args = parser.parse_args()

    np.random.seed(args.seed)
    out = Path(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)

    groups = [g.strip() for g in args.groups.split(",")]
    n_per_group = args.n_samples // len(groups)
    n_otus = min(args.n_otus, len(TAXA_LIBRARY))
    selected_taxa = TAXA_LIBRARY[:n_otus]

    all_counts = []
    all_ids    = []
    all_groups = []

    for g in groups:
        mods = GROUP_MODIFIERS.get(g, DEFAULT_MODIFIER)
        baseline = np.array([
            row[4] * mods.get(row[3], 1.0) for row in selected_taxa
        ], dtype=float)
        counts = simulate_group(baseline, n_per_group)
        prefix = g[0].upper()
        ids = [f"{prefix}{i+1:02d}" for i in range(n_per_group)]
        all_counts.append(counts)
        all_ids.extend(ids)
        all_groups.extend([g] * n_per_group)

    otu_array = np.vstack(all_counts)
    otu_ids   = [row[0] for row in selected_taxa]

    otu_table = pd.DataFrame(otu_array, index=all_ids, columns=otu_ids)
    otu_table.index.name = "sample_id"

    taxonomy = pd.DataFrame(
        [(row[0], row[1], row[2], row[3]) for row in selected_taxa],
        columns=["OTU_ID","genus","family","phylum"]
    ).set_index("OTU_ID")

    metadata = pd.DataFrame({
        "sample_id": all_ids,
        "group":     all_groups,
        "subject":   [f"S{i+1:02d}" for i in range(len(all_ids))],
        "age":       np.random.randint(22, 55, len(all_ids)),
        "sex":       np.random.choice(["M","F"], len(all_ids)),
    }).set_index("sample_id")

    otu_table.to_csv(out / "otu_table.csv")
    taxonomy.to_csv(out / "taxonomy.csv")
    metadata.to_csv(out / "metadata.csv")

    print(f"✓ Generated synthetic data in {out}/")
    print(f"  {len(all_ids)} samples × {n_otus} OTUs")
    for g in groups:
        print(f"  {g}: n={all_groups.count(g)}")


if __name__ == "__main__":
    main()
