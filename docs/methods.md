# Methods notes: plastizyme screen

- Input: 24 field MAGs (YJ 13, JJ 11), ≥50 % completeness and ≤10 % contamination (CheckM v1.0.18).
- Gene prediction: Prodigal v2.6.3 `-p single` (fallback `-p meta`); protein IDs prefixed `SAMPLE__MAG__`.
- References: PAZy and PlasticDB PLA/PBAT enzymes (files in `references/`); merged per polymer,
  sequences <100 aa removed, CD-HIT v4.8.1 `-c 0.95 -n 5`.
- Search: DIAMOND v2.1.10 `blastp --more-sensitive --evalue 1e-5 --id 30 --query-cover 70`.
- Counting: a predicted gene is a putative plastizyme if it has ≥1 hit to either reference set; genes hitting
  both sets are counted once in "any hit". Normalised per MAG and per 1,000 predicted genes.
- Same settings as the plastizyme module of darsha-crypto/bioplastic-soil-multiomics-pipeline
  (`modules/plastizyme_qc.nf`, `modules/diamond.nf`).
