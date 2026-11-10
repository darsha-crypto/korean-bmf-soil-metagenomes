# Genome-resolved metagenomic dataset of Korean agricultural soils exposed to PLA–PBAT biodegradable mulch film

Analysis code and derived tables accompanying the Data Descriptor:

> Prabhaharan D., Lee Y., Yoon G., Chung J., Ahn J.-H., Park H. & Sang B.-I.
> *Genome-resolved metagenomic dataset of Korean agricultural soils exposed to PLA–PBAT biodegradable mulch film.* (in preparation)

The dataset comprises seven shotgun metagenomes and 177 metagenome-assembled genomes (MAGs) from two South Korean
agricultural soils with contrasting biodegradable mulch film (BMF) histories, and from an 80-day laboratory
incubation of one of them with powdered PLA–PBAT BMF. This repository holds the plastizyme screen of the field MAGs
(Table S5 of the manuscript) and the scripts to reproduce it.

## Samples

| Sample | Description | BMF | MAGs |
|---|---|---|---|
| YJ | Field soil, Yeoju (no BMF-use history) | – | 13 |
| JJ | Field soil, Jeonju (multi-year BMF use) | – | 11 |
| NC | Natural JJ soil, incubation control | – | 15 |
| NP1, NP2 | Natural JJ soil + powdered BMF | + | 23, 37 |
| SC | Standard soil mixture (modified ISO 17556), control | – | 34 |
| SP | Standard soil mixture + powdered BMF | + | 44 |

Incubation: 80 days, 30 °C, in the dark, 500 g soil in 2.4-L reactors, 40–50 % of water-holding capacity,
soil-to-BMF ratio 250:1 (w/w).

## Data availability

| Data | Repository | Accession |
|---|---|---|
| Raw reads | NCBI SRA | BioProject PRJNA1134933; SRR29792059–SRR29792065 |
| Per-sample assemblies | GenBank WGS | JBJSWA000000000–JBJSWG000000000 |
| MAGs, annotations, DRAM and gene tables | figshare | DOI to be added |

MAG sequences are distributed through figshare and NCBI, not through this repository.

## Plastizyme screen of the field MAGs

Putative PLA- and PBAT-degrading enzymes were searched in the 24 MAGs from the two field soils (YJ, JJ).

**Method** (`scripts/run_plastizyme_screen.sh`)

1. **Gene prediction:** Prodigal v2.6.3, single-genome mode (metagenome mode as fallback).
2. **Reference enzymes:** experimentally characterised PLA and PBAT depolymerases from PAZy and PlasticDB
   (`references/`), merged per polymer, filtered to ≥100 aa and de-duplicated with CD-HIT v4.8.1 at 95 % identity
   (PLA 69 → 46 sequences; PBAT 38 → 27 sequences; `results/qc/`).
3. **Search:** DIAMOND v2.1.10 blastp `--more-sensitive`, E ≤ 1e-5, identity ≥ 30 %, query coverage ≥ 70 %.
4. **Summary:** `scripts/summarize_plastizyme_hits.py`; a gene is counted once even if it hits both reference sets.

**Result** (`results/plastizyme_summary_per_sample.tsv`)

| Sample | MAGs | MAGs with ≥1 hit | Predicted genes | Genes with PLA hits | Genes with PBAT hits | Genes with any hit | Per MAG | Per 1,000 genes |
|---|---|---|---|---|---|---|---|---|
| YJ (no BMF history) | 13 | 13 | 42,110 | 41 | 21 | 50 | 3.85 | 1.19 |
| JJ (multi-year BMF use) | 11 | 10 | 43,691 | 64 | 31 | 72 | 6.55 | 1.65 |

MAGs from the soil with a BMF-use history (JJ) encoded more putative PLA- and PBAT-degrading enzymes than those
from the soil without (YJ), whether normalised per MAG or per 1,000 predicted genes. Hits are homology-based
(≥30 % identity, mostly carboxylesterases, lipases and other α/β-hydrolases) and should be read as putative
plastizymes; MAGs capture only part of each community.

## Repository layout

```
references/   PAZy and PlasticDB PLA/PBAT reference proteins (input)
scripts/      run_plastizyme_screen.sh, summarize_plastizyme_hits.py
notebooks/    Google Colab version of the screen
results/
  qc/                               reference QC reports and de-duplicated reference sets
  diamond/                          raw DIAMOND hit tables (PLA, PBAT)
  plastizyme_summary_per_sample.tsv per-sample counts and normalised values
  TableS5_plastizyme_per_MAG.tsv    per-MAG counts (Table S5)
  plastizyme_hits_per_gene.tsv      best reference hit for every positive gene
docs/         methods notes
```

## Reproducing the screen

Locally (Linux or macOS):

```bash
conda env create -f environment.yml
conda activate bmf-plastizyme
# MAG FASTA files from figshare, one folder per sample: bins/YJ/*.fasta, bins/JJ/*.fasta
bash scripts/run_plastizyme_screen.sh bins references screen_out 4
python scripts/summarize_plastizyme_hits.py bins screen_out results
```

In Google Colab: open `notebooks/plastizyme_screen_colab.ipynb` with **File → Open notebook → GitHub**.

Assembly, binning and annotation (FastQC, Trimmomatic, MEGAHIT, metaSPAdes, QUAST, MetaBAT2, MaxBin2, CONCOCT,
DAS Tool, CheckM, GTDB-Tk, RASTtk, SpeciesTree, dbCAN2, DRAM) were run in KBase; software versions and
non-default parameters are given in the manuscript Methods.

## Citation

Please cite the Data Descriptor (see `CITATION.cff`) and the reference databases used:
PAZy (Buchholz et al., *Proteins* 2022) and PlasticDB (Gambarini et al., *Database* 2022).

## License

Code: MIT (see `LICENSE`). Reference sequences in `references/` remain subject to the terms of PAZy and PlasticDB.
