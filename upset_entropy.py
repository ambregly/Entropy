#!/usr/bin/env python3
"""Compare les ensembles d'identifiants de fusion entre les filtres d'entropie
k3, k5, k8 et k10, et produit deux UpSet plots : un pour les kmers, un pour les
contigs.

Fichiers attendus dans le dossier `complexity` (voir --input) :
    kmers-entropyk3.fa   kmers-entropyk5.fa   kmers-entropyk8.fa   kmers-entropyk10.fa
    contigs-entropyk3.fa contigs-entropyk5.fa contigs-entropyk8.fa contigs-entropyk10.fa

Format des en-tetes FASTA :
    >ELOVL5_chr6_53213614_53213640_PTP4A1_chr6_64286340_64286365_2744.kmer1
    >OAZ1_chr19_2269717_2269743_ZBTB7A_chr19_4055220_4055245_2564.contig_1

La cle de comparaison est l'IDENTIFIANT UNIQUE de fusion, c.-a-d. l'en-tete
prive du suffixe `.kmerN` / `.contig_N`. Une fusion est consideree presente pour
un k donne des qu'au moins un kmer (resp. contig) la concerne dans le fichier
correspondant.

Sorties (dans --outdir) :
    upset_kmers.png / .pdf        UpSet plot des fusions presentes par k (kmers)
    upset_contigs.png / .pdf      idem pour les contigs
    membership_kmers.csv          tableau presence/absence (fusion x k)
    membership_contigs.csv        idem pour les contigs

Dependances : pip install upsetplot pandas matplotlib

Exemple :
    python upset_entropy.py --input ~/complexity --outdir ~/complexity/upset
"""

import argparse
import os
import re
import sys
from glob import glob

# Ordre canonique des valeurs de k a afficher dans les UpSet plots.
K_ORDER = [3, 5, 8, 10]

# Capture la valeur de k dans un nom de fichier du type "kmers-entropyk10.fa".
K_IN_FILENAME = re.compile(r"entropy[_-]?k(\d+)", re.IGNORECASE)

# Supprime le suffixe ".kmerN" ou ".contig_N" (avec ou sans underscore) en fin
# d'en-tete pour ne garder que l'identifiant unique de fusion.
SUFFIX = re.compile(r"\.(?:kmer|contig)s?_?\d+\s*$", re.IGNORECASE)


def extract_id(header):
    """Retourne l'identifiant unique de fusion a partir d'une ligne d'en-tete."""
    header = header[1:].strip()          # retire le '>' initial
    header = header.split()[0]           # ignore un eventuel commentaire apres espace
    return SUFFIX.sub("", header)


def read_ids(path):
    """Lit un FASTA et renvoie l'ensemble des identifiants uniques de fusion."""
    ids = set()
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                ids.add(extract_id(line))
    return ids


def discover(input_dir, prefix):
    """Associe chaque valeur de k a son fichier `{prefix}-entropyk*.fa`.

    Renvoie un dict {k: chemin} ordonne selon K_ORDER pour les k presents.
    """
    found = {}
    for path in glob(os.path.join(input_dir, f"{prefix}-entropy*.fa")):
        m = K_IN_FILENAME.search(os.path.basename(path))
        if m:
            found[int(m.group(1))] = path
    return {k: found[k] for k in K_ORDER if k in found}


def build_membership(files_by_k):
    """Construit un DataFrame booleen (fusion x k) de presence/absence."""
    import pandas as pd

    ids_by_k = {k: read_ids(p) for k, p in files_by_k.items()}
    all_ids = sorted(set().union(*ids_by_k.values())) if ids_by_k else []
    columns = [f"k{k}" for k in files_by_k]
    data = {
        f"k{k}": [fid in ids_by_k[k] for fid in all_ids]
        for k in files_by_k
    }
    df = pd.DataFrame(data, index=all_ids, columns=columns)
    df.index.name = "fusion_id"
    return df


def make_upset(df, title, out_base):
    """Genere et sauvegarde un UpSet plot a partir d'un DataFrame booleen."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from upsetplot import UpSet, from_indicators

    # from_indicators construit la serie multi-index attendue par UpSet a partir
    # des colonnes booleennes (une par k).
    data = from_indicators(list(df.columns), data=df)

    # NB: on n'utilise PAS show_counts=True : l'annotation interne d'upsetplot
    # 0.9.0 est incompatible avec numpy >= 2 ("only 0-dimensional arrays can be
    # converted to Python scalars"). On ajoute donc les compteurs nous-memes via
    # bar_label, ce qui fonctionne quelle que soit la version de numpy.
    upset = UpSet(
        data,
        show_counts=False,
        sort_by="cardinality",
        sort_categories_by="-input",  # garde l'ordre k3, k5, k8, k10
        min_subset_size=1,
    )
    axes = upset.plot()

    # Etiquette chaque barre avec son effectif (intersections + totaux par k).
    for key in ("intersections", "totals"):
        ax = axes.get(key)
        if ax is None:
            continue
        for container in ax.containers:
            ax.bar_label(container, fmt="%d", padding=2, fontsize=8)

    plt.suptitle(title)
    for ext in ("png", "pdf"):
        out = f"{out_base}.{ext}"
        plt.savefig(out, dpi=200, bbox_inches="tight")
        print(f"  -> {out}")
    plt.close("all")


def process(kind, prefix, input_dir, outdir):
    """Traite un type ('kmers' ou 'contigs') : fichiers -> membership -> UpSet."""
    files_by_k = discover(input_dir, prefix)
    if not files_by_k:
        print(f"[{kind}] aucun fichier '{prefix}-entropyk*.fa' dans {input_dir} "
              f"- ignore.")
        return
    print(f"[{kind}] fichiers detectes :")
    for k, p in files_by_k.items():
        print(f"    k{k}: {os.path.basename(p)}")

    df = build_membership(files_by_k)
    print(f"[{kind}] {len(df)} identifiants uniques au total.")
    for col in df.columns:
        print(f"    {col}: {int(df[col].sum())} fusions")

    csv_path = os.path.join(outdir, f"membership_{kind}.csv")
    df.astype(int).to_csv(csv_path)
    print(f"  -> {csv_path}")

    make_upset(df, f"UpSet des {kind} par filtre d'entropie (k)",
               os.path.join(outdir, f"upset_{kind}"))


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--input", default=os.path.expanduser("~/complexity"),
                        help="dossier contenant les fichiers .fa "
                             "(defaut: ~/complexity)")
    parser.add_argument("--outdir", default=None,
                        help="dossier de sortie (defaut: <input>/upset)")
    args = parser.parse_args(argv)

    input_dir = os.path.expanduser(args.input)
    if not os.path.isdir(input_dir):
        sys.exit(f"Erreur: dossier introuvable: {input_dir}")
    outdir = args.outdir or os.path.join(input_dir, "upset")
    os.makedirs(outdir, exist_ok=True)

    process("kmers", "kmers", input_dir, outdir)
    process("contigs", "contigs", input_dir, outdir)
    print("Termine.")


if __name__ == "__main__":
    main()
