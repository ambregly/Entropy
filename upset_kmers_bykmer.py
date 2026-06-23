#!/usr/bin/env python3
"""UpSet plot comparant les KMERS INDIVIDUELS (et non les fusions) entre les
filtres d'entropie k3, k5, k8 et k10.

Contrairement a upset_entropy.py qui compare les identifiants de fusion, ici
chaque kmer est identifie par son EN-TETE COMPLET `fusion.kmerN` (ex:
ELOVL5_chr6_..._2744.kmer1). Deux kmers de fusions differentes restent donc
distincts meme s'ils partagent la meme sequence.

Fichiers attendus dans le dossier `complexity` (voir --input) :
    kmers-entropyk{k}.fa  (ou kmersentropyk{k}.fa)

Sorties (dans --outdir) :
    upset_kmers_bykmer.png / .pdf   UpSet plot des kmers par k
    membership_kmers_bykmer.csv     tableau presence/absence (kmer x k)

Dependances : voir requirements.txt

Exemple :
    python upset_kmers_bykmer.py --input ~/complexity --outdir ~/complexity/upset
"""

import argparse
import os
import sys

# Reutilise la detection de fichiers et le trace UpSet du script principal.
from upset_entropy import K_ORDER, discover, make_upset


def read_kmer_headers(path):
    """Renvoie l'ensemble des en-tetes complets de kmers (`fusion.kmerN`).

    On retire le '>' initial et l'eventuel commentaire apres l'espace
    (ex: ` ct:1`), en gardant le 1er champ = l'identifiant complet du kmer.
    """
    keys = set()
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                keys.add(line[1:].split()[0])
    return keys


def build_membership(files_by_k):
    """DataFrame booleen (kmer x k) de presence/absence des kmers individuels."""
    import pandas as pd

    keys_by_k = {k: read_kmer_headers(p) for k, p in files_by_k.items()}
    all_keys = sorted(set().union(*keys_by_k.values())) if keys_by_k else []
    data = {
        f"k{k}": [key in keys_by_k[k] for key in all_keys]
        for k in files_by_k
    }
    df = pd.DataFrame(data, index=all_keys, columns=[f"k{k}" for k in files_by_k])
    df.index.name = "kmer_id"
    return df


def main(argv=None):
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--input", default=os.path.expanduser("~/complexity"),
                        help="dossier contenant les fichiers .fa (defaut: ~/complexity)")
    parser.add_argument("--outdir", default=None,
                        help="dossier de sortie (defaut: <input>/upset)")
    args = parser.parse_args(argv)

    input_dir = os.path.expanduser(args.input)
    if not os.path.isdir(input_dir):
        sys.exit(f"Erreur: dossier introuvable: {input_dir}")
    outdir = args.outdir or os.path.join(input_dir, "upset")
    os.makedirs(outdir, exist_ok=True)

    files_by_k = discover(input_dir, "kmers")
    if not files_by_k:
        sys.exit(f"Aucun fichier 'kmers*entropyk*.fa' dans {input_dir}.")
    print("Fichiers detectes :")
    for k, p in files_by_k.items():
        print(f"    k{k}: {os.path.basename(p)}")

    df = build_membership(files_by_k)
    print(f"{len(df)} kmers uniques au total.")
    for col in df.columns:
        print(f"    {col}: {int(df[col].sum())} kmers")

    csv_path = os.path.join(outdir, "membership_kmers_bykmer.csv")
    df.astype(int).to_csv(csv_path)
    print(f"  -> {csv_path}")

    if len(df.columns) < 2:
        print(f"Un seul k present ({', '.join(df.columns)}) : pas d'UpSet plot "
              f"(comparaison impossible avec un seul ensemble). CSV produit.")
        return

    make_upset(df, "UpSet des kmers (par kmer) par filtre d'entropie (k)",
               os.path.join(outdir, "upset_kmers_bykmer"))
    print("Termine.")


if __name__ == "__main__":
    main()
