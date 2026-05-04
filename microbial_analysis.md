---
name: microbiome-analysis
description: >
  Automatically generate a complete microbiome analysis pipeline in Python from an OTU/feature
  table (or raw count matrix). Use this skill whenever the user mentions microbiome, 16S, amplicon,
  metagenomics, OTU table, taxa, alpha diversity, beta diversity, Shannon, Bray-Curtis, QIIME,
  or asks to analyze bacterial community data. Also triggers when the user uploads a CSV that
  looks like a samples × species count matrix, or asks to compare microbial communities between
  groups. Always prefer this skill over writing ad hoc microbiome code — it produces consistent,
  well-documented, publication-quality outputs every time.
---

# Microbiome Analysis Skill

Generates a complete, documented microbiome analysis pipeline. Handles both **synthetic data**
(for learning/demos) and **real data** (user-uploaded OTU tables, QIIME2 exports, Bracken outputs).

## When to use this skill

- User asks to analyze a microbiome dataset
- User uploads an OTU table, feature table, or taxa count matrix
- User wants to compare microbial diversity between groups (healthy vs disease, pre vs post treatment, etc.)
- User wants to learn how 16S or metagenomics analysis works
- User has Kraken2/Bracken output and wants diversity metrics + figures

---

## Step 0: Understand the inputs

Before writing any code, determine:

1. **Data source** — does the user have real data, or do they want a synthetic demo?
   - Real data: ask for OTU table, taxonomy file, metadata file (or QIIME2 `.qza`)
   - No data: generate synthetic data using `scripts/generate_synthetic.py` (see below)

2. **Input format** — read the file to determine:
   - Rows = samples or OTUs? (rows should be samples for analysis; transpose if needed)
   - Is there a taxonomy file? (maps OTU IDs → genus, family, phylum)
   - Is there a metadata file? (maps sample IDs → group labels, covariates)

3. **Comparison groups** — what groups exist in the metadata? (e.g. Healthy vs Disease)

4. **Analysis depth** — what does the user want?
   - Quick: taxonomy barplot only
   - Standard: taxonomy + alpha diversity + beta diversity PCoA ← **default**
   - Extended: + statistical testing, rarefaction curves, differential abundance

Read `references/input_formats.md` if you're unsure how to parse the user's file format.

---

## Step 1: Set up the environment

Always run this first:

```bash
pip install pandas numpy scipy matplotlib seaborn scikit-learn --break-system-packages -q
```

No conda environments needed — all standard PyPI packages.

---

## Step 2: Generate or load data

**If generating synthetic data** (no user files):

```bash
python "C:\Users\ninag\.claude\commands\scripts_for_microbial_analysis\generate_synthetic.py" \
  --n-samples 20 \
  --n-otus 20 \
  --groups "Healthy,Antibiotic" \
  --output-dir .
```

This creates: `otu_table.csv`, `taxonomy.csv`, `metadata.csv`

**If user provides real data**, load and validate:

```python
import pandas as pd

otu = pd.read_csv("otu_table.csv", index_col=0)
# Ensure rows = samples, columns = OTUs
if otu.shape[0] > otu.shape[1]:
    print("WARNING: More rows than columns — check if transposition is needed")
```

Read `references/input_formats.md` for QIIME2 `.qza` and Bracken formats.

---

## Step 3: Run the analysis

Use `scripts/run_analysis.py` as the base. It accepts arguments so you don't need to
rewrite it for each user — just pass the right flags:

```bash
python "C:\Users\ninag\.claude\commands\scripts_for_microbial_analysis\run_analysis.py" \
  --otu otu_table.csv \
  --taxonomy taxonomy.csv \
  --metadata metadata.csv \
  --group-col group \
  --output-dir results/
```

### What the script produces

| Output file | Contents |
|---|---|
| `results/microbiome_analysis.png` | 3-panel figure: taxonomy barplot, alpha diversity boxplot, beta diversity PCoA |
| `results/alpha_diversity.csv` | Shannon entropy, Simpson index, richness per sample |
| `results/beta_diversity_matrix.csv` | Pairwise Bray-Curtis distance matrix |
| `results/phylum_abundance.csv` | Phylum-level relative abundance per sample |

---

## Step 4: Interpret and explain the results

After generating figures, always explain the three panels to the user:

**Taxonomy barplot:**
- What phyla dominate in each group?
- Look for: Proteobacteria bloom (dysbiosis), loss of Firmicutes/Bacteroidetes (antibiotic effect)
- Key ratio: Firmicutes:Bacteroidetes (F:B ratio) — disrupted in many disease states

**Alpha diversity:**
- Report the Mann-Whitney U p-value
- Higher Shannon = more diverse, more even community
- Lower Shannon = dominated by one or few taxa (often a sign of dysbiosis)
- Richness = raw count of taxa present (ignores evenness)

**Beta diversity PCoA:**
- Report within-group vs between-group Bray-Curtis distances
- Clear separation = groups have distinct community composition
- Overlap = groups are similar

---

## Step 5: Add a README and push to GitHub (if requested)

Use `scripts/generate_readme.py` to auto-generate a README with:
- Dataset description
- Results figure embedded
- Methods summary
- Quickstart instructions

```bash
python "C:\Users\ninag\.claude\commands\scripts_for_microbial_analysis\generate_readme.py" \
  --title "Gut microbiome analysis" \
  --dataset-desc "Healthy vs antibiotic-treated subjects" \
  --figure results/microbiome_analysis.png \
  --output README.md
```

---

## Decision table

| User says | What to do |
|---|---|
| "Analyze my OTU table" | Load their file → Step 3 → Step 4 |
| "Show me how microbiome analysis works" | Synthetic data → all steps with explanations |
| "Just the figure, no code explanation" | Run script, present figure, brief interpretation only |
| "I have QIIME2 output" | Read `references/input_formats.md` first |
| "I have Bracken output" | Read `references/input_formats.md` → map to OTU table format |
| "Add statistical testing" | Read `references/stats.md` for PERMANOVA + differential abundance |

---

## File map

```
microbiome-analysis-skill/
├── SKILL.md                        ← you are here
├── scripts/
│   ├── generate_synthetic.py       ← generates demo data
│   ├── run_analysis.py             ← main analysis pipeline
│   └── generate_readme.py          ← auto-generates GitHub README
└── references/
    ├── input_formats.md            ← how to parse QIIME2, Bracken, biom formats
    └── stats.md                    ← PERMANOVA, differential abundance, rarefaction
```

Read reference files only when needed (e.g. user has QIIME2 data → read `input_formats.md`).
Default runs don't require them.
