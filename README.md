# Gut Microbiome Analysis — Healthy vs. Antibiotic-Treated

> **Prototype:** This is an early-stage project. I am actively learning microbial analysis and plan to expand and improve this pipeline as my knowledge grows.

## Overview

This project runs a basic microbiome diversity analysis comparing two groups: **Healthy** and **Antibiotic-treated** samples. It produces a three-panel figure covering community composition, alpha diversity, and beta diversity.

The current example uses synthetic data generated for learning and development purposes but this claude skill can be used on real data.

## Output

The analysis generates:

| File | Description |
|---|---|
| `results/microbiome_analysis.png` | Taxonomy barplot, alpha diversity boxplot, beta diversity PCoA |
| `results/alpha_diversity.csv` | Shannon entropy, Simpson index, and richness per sample |
| `results/beta_diversity_matrix.csv` | Pairwise Bray-Curtis distance matrix |
| `results/phylum_abundance.csv` | Phylum-level relative abundance per sample |

## How to Run

**Install dependencies:**
```bash
pip install pandas numpy scipy matplotlib seaborn scikit-learn
```

**Generate synthetic data:**
```bash
python scripts_for_microbial_analysis/generate_synthetic.py \
  --n-samples 20 \
  --n-otus 20 \
  --groups "Healthy,Antibiotic" \
  --output-dir .
```

**Run the analysis:**
```bash
python scripts_for_microbial_analysis/run_analysis.py \
  --otu otu_table.csv \
  --taxonomy taxonomy.csv \
  --metadata metadata.csv \
  --group-col group \
  --output-dir results/
```

## What the Results Show

- **Taxonomy barplot** — relative abundance of bacterial phyla per sample, grouped by condition
- **Alpha diversity** — within-sample diversity (Shannon entropy); higher = more even community
- **Beta diversity PCoA** — between-sample diversity using Bray-Curtis distances; separated clusters indicate distinct community compositions

## Planned Improvements

- Replace synthetic data with real 16S rRNA sequencing data
- Add statistical testing (PERMANOVA) for beta diversity
- Add rarefaction curves to assess sequencing depth
- Add differential abundance analysis
- Improve phylum-level visualizations

## Attribution

Pipeline scripts adapted from the [microbiome-analysis Claude Code skill](https://github.com/anthropics/claude-code).
