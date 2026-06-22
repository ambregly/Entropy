#!/usr/bin/env python3
"""Construit un tableau du NOMBRE DE KMERS associes a chaque identifiant unique
de fusion, pour chaque filtre d'entropie (k3, k5, k8, k10).

Fichiers attendus dans le dossier `complexity` (voir --input) :
    kmers-entropyk{k}.fa  (ou kmersentropyk{k}.fa)

Chaque en-tete `>FUSION_ID.kmerN ...` compte pour 1 kmer de la fusion FUSION_ID
(l'identifiant unique = en-tete prive du suffixe `.kmerN`). Le tableau donne,
par fusion et par k, le nombre de lignes `.kmerN` presentes (0 si absente).

Sorties (dans --outdir) :
    kmer_counts.csv   colonnes: fusion_id, k3, k5, k8, k10 (nombre de kmers)

Exemple :
    python kmer_counts.py --input ~/complexity --outdir ~/complexity
"""

import argparse
import os
import sys
from collections import Counter

# Reutilise la logique d'extraction d'identifiant et de detection de fichiers.
from upset_entropy import K_ORDER, discover, extract_id


def count_ids(path):
    """Renvoie un Counter {fusion_id: nombre de kmers} pour un fichier FASTA."""
    counts = Counter()
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                counts[extract_id(line)] += 1
    return counts


def build_table(files_by_k):
    """DataFrame (fusion_id x k) du nombre de kmers, 0 si la fusion est absente."""
    import pandas as pd

    counts_by_k = {k: count_ids(p) for k, p in files_by_k.items()}
    all_ids = sorted(set().union(*(c.keys() for c in counts_by_k.values())))
    data = {
        f"k{k}": [counts_by_k[k].get(fid, 0) for fid in all_ids]
        for k in files_by_k
    }
    df = pd.DataFrame(data, index=all_ids,
                      columns=[f"k{k}" for k in files_by_k])
    df.index.name = "fusion_id"
    return df


def main(argv=None):
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--input", default=os.path.expanduser("~/complexity"),
                        help="dossier contenant les fichiers .fa (defaut: ~/complexity)")
    parser.add_argument("--outdir", default=None,
                        help="dossier de sortie (defaut: = --input)")
    args = parser.parse_args(argv)

    input_dir = os.path.expanduser(args.input)
    if not os.path.isdir(input_dir):
        sys.exit(f"Erreur: dossier introuvable: {input_dir}")
    outdir = args.outdir or input_dir
    os.makedirs(outdir, exist_ok=True)

    files_by_k = discover(input_dir, "kmers")
    if not files_by_k:
        sys.exit(f"Aucun fichier 'kmers*entropyk*.fa' dans {input_dir}.")
    print("Fichiers detectes :")
    for k, p in files_by_k.items():
        print(f"    k{k}: {os.path.basename(p)}")

    df = build_table(files_by_k)
    out = os.path.join(outdir, "kmer_counts.csv")
    df.to_csv(out)
    print(f"{len(df)} identifiants uniques.")
    for col in df.columns:
        present = int((df[col] > 0).sum())
        total = int(df[col].sum())
        print(f"    {col}: {present} fusions, {total} kmers au total "
              f"(moyenne {total / present:.1f} kmers/fusion)" if present else
              f"    {col}: 0 fusion")
    print(f"  -> {out}")


if __name__ == "__main__":
    main()
