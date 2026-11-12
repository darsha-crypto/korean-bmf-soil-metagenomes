#!/usr/bin/env python3
"""Summarise DIAMOND plastizyme hits per gene, per MAG and per sample.

Usage:
    summarize_plastizyme_hits.py <bins_root> <screen_outdir> <results_dir>

<bins_root>/<SAMPLE>/*.fa|*.fasta   MAG nucleotide FASTA, one sub-folder per sample
<screen_outdir>                     output folder of run_plastizyme_screen.sh
<results_dir>                       where the three summary tables are written

A gene counts as a putative plastizyme if it has at least one DIAMOND hit
(E <= 1e-5, identity >= 30 %, query coverage >= 70 %) against the PLA or PBAT
reference set. Per-sample values are normalised per MAG and per 1,000 genes.
"""
import csv
import glob
import os
import sys
from collections import defaultdict

DBS = ("PLA", "PBAT")


def fasta_length(path):
    with open(path) as fh:
        return sum(len(line.strip()) for line in fh if not line.startswith(">"))


def count_records(path):
    with open(path) as fh:
        return sum(1 for line in fh if line.startswith(">"))


def main(bins_root, screen_out, results):
    os.makedirs(results, exist_ok=True)

    # MAG inventory: id -> (sample, MAG, size bp, number of predicted genes)
    mags = {}
    for sample_dir in sorted(glob.glob(os.path.join(bins_root, "*", ""))):
        sample = os.path.basename(os.path.dirname(sample_dir))
        for fasta in sorted(glob.glob(sample_dir + "*.fa") + glob.glob(sample_dir + "*.fasta")):
            mag = os.path.splitext(os.path.basename(fasta))[0]
            mag_id = f"{sample}__{mag}"
            faa = os.path.join(screen_out, "prodigal", f"{mag_id}.faa")
            mags[mag_id] = (sample, mag, fasta_length(fasta), count_records(faa))

    # Reference descriptions (first header token -> full header)
    ref_desc = {}
    for db in DBS:
        with open(os.path.join(screen_out, "qc", f"{db}.clean.faa")) as fh:
            for line in fh:
                if line.startswith(">"):
                    tok = line[1:].split()[0]
                    ref_desc[tok] = line[1:].strip()

    # DIAMOND hits: best hit (highest bit score) per gene and database
    best = {}
    for db in DBS:
        with open(os.path.join(screen_out, "diamond", f"{db}.diamond.tsv")) as fh:
            for row in csv.reader(fh, delimiter="\t"):
                gene, ref = row[0], row[1]
                pident, evalue, bits, qcov = float(row[2]), row[10], float(row[11]), row[12]
                key = (gene, db)
                if key not in best or bits > best[key][3]:
                    best[key] = (ref, pident, evalue, bits, qcov)

    # Per-gene table
    genes = defaultdict(dict)
    for (gene, db), hit in best.items():
        genes[gene][db] = hit
    with open(os.path.join(results, "plastizyme_hits_per_gene.tsv"), "w", newline="") as out:
        w = csv.writer(out, delimiter="\t")
        w.writerow(["sample", "MAG", "gene", "database", "best_reference", "identity_pct",
                    "evalue", "bitscore", "query_coverage_pct", "reference_description"])
        for gene in sorted(genes):
            sample, mag = gene.split("__")[:2]
            for db in DBS:
                if db in genes[gene]:
                    ref, pid, ev, bits, qcov = genes[gene][db]
                    w.writerow([sample, mag, gene, db, ref, f"{pid:.1f}", ev, f"{bits:.0f}", qcov,
                                ref_desc.get(ref, "")])

    # Per-MAG table (Table S2)
    per_mag = defaultdict(lambda: {db: set() for db in DBS})
    for (gene, db) in best:
        per_mag["__".join(gene.split("__")[:2])][db].add(gene)
    mag_rows = []
    for mag_id, (sample, mag, bp, n_genes) in sorted(mags.items()):
        pla, pbat = per_mag[mag_id]["PLA"], per_mag[mag_id]["PBAT"]
        mag_rows.append([sample, mag, round(bp / 1e6, 2), n_genes, len(pla), len(pbat), len(pla | pbat)])
    with open(os.path.join(results, "TableS2_plastizyme_per_MAG.tsv"), "w", newline="") as out:
        w = csv.writer(out, delimiter="\t")
        w.writerow(["sample", "MAG", "size_Mb", "predicted_genes", "PLA_hit_genes",
                    "PBAT_hit_genes", "any_hit_genes"])
        w.writerows(mag_rows)

    # Per-sample summary
    with open(os.path.join(results, "plastizyme_summary_per_sample.tsv"), "w", newline="") as out:
        w = csv.writer(out, delimiter="\t")
        header = ["sample", "MAGs", "MAGs_with_hit", "predicted_genes", "MAG_size_Mb",
                  "PLA_hit_genes", "PBAT_hit_genes", "any_hit_genes",
                  "hit_genes_per_MAG", "hit_genes_per_1000_genes", "hit_genes_per_Mb"]
        w.writerow(header)
        print("\t".join(header))
        for sample in sorted({r[0] for r in mag_rows}):
            rows = [r for r in mag_rows if r[0] == sample]
            n_genes = sum(r[3] for r in rows)
            size = sum(r[2] for r in rows)
            anyhit = sum(r[6] for r in rows)
            line = [sample, len(rows), sum(1 for r in rows if r[6] > 0), n_genes, round(size, 1),
                    sum(r[4] for r in rows), sum(r[5] for r in rows), anyhit,
                    round(anyhit / len(rows), 2), round(anyhit / n_genes * 1000, 2),
                    round(anyhit / size, 2)]
            w.writerow(line)
            print("\t".join(map(str, line)))


if __name__ == "__main__":
    if len(sys.argv) != 4:
        sys.exit(__doc__)
    main(*sys.argv[1:])
