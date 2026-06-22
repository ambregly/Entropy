#!/usr/bin/env python3
"""Pour chaque filtre d'entropie (k3, k5, k8, k10), trace un diagramme de Venn
comparant les identifiants uniques de fusion presents dans les kmers vs dans les
contigs.

Fichiers attendus dans le dossier `complexity` (voir --input) :
    kmers-entropyk{k}.fa     et     contigs-entropyk{k}.fa

La cle de comparaison est l'IDENTIFIANT UNIQUE de fusion (en-tete prive du
suffixe `.kmerN` / `.contig_N`), exactement comme dans upset_entropy.py.

Sorties (dans --outdir) :
    venn_k{k}.png / .pdf      un diagramme de Venn kmers vs contigs par k
    venn_all.png / .pdf       panneau recapitulatif (tous les k cote a cote)
    venn_counts.csv           effectifs (kmers_seuls, contigs_seuls, intersection) par k

Dependances : pip install matplotlib-venn (+ celles de requirements.txt)

Exemple :
    python venn_entropy.py --input ~/complexity --outdir ~/complexity/venn
"""

import argparse
import math
import os
import sys

# Reutilise la logique de lecture/extraction d'identifiants du script UpSet.
from upset_entropy import K_ORDER, discover, read_ids


def venn_pair(ax, kmer_ids, contig_ids, k, show_title=True):
    """Trace un Venn 2 ensembles (kmers vs contigs) sur l'axe donne."""
    from matplotlib_venn import venn2

    v = venn2(
        [kmer_ids, contig_ids],
        set_labels=("kmers", "contigs"),
        ax=ax,
    )
    # Quand un ensemble est presque inclus dans l'autre, les deux etiquettes
    # d'ensemble se superposent sous les cercles ("kmerscontigs"). On les
    # repositionne en haut a gauche / en haut a droite pour eviter la collision.
    labels = v.set_labels
    if labels and labels[0] is not None:
        labels[0].set_position((-0.6, 0.65))
        labels[0].set_horizontalalignment("center")
    if labels and labels[1] is not None:
        labels[1].set_position((0.6, 0.65))
        labels[1].set_horizontalalignment("center")
    if show_title:
        ax.set_title(f"k{k}")
    return v


def process(input_dir, outdir):
    """Pour chaque k disponible, lit les deux fichiers et trace un Venn."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    kmer_files = discover(input_dir, "kmers")
    contig_files = discover(input_dir, "contigs")
    ks = [k for k in K_ORDER if k in kmer_files and k in contig_files]
    if not ks:
        sys.exit("Aucun couple kmers/contigs trouve "
                 f"(kmers: {sorted(kmer_files)}, contigs: {sorted(contig_files)}).")

    # Lecture des ensembles d'identifiants par k.
    ids = {}
    rows = []
    for k in ks:
        kmer_ids = read_ids(kmer_files[k])
        contig_ids = read_ids(contig_files[k])
        ids[k] = (kmer_ids, contig_ids)
        inter = len(kmer_ids & contig_ids)
        rows.append((k, len(kmer_ids - contig_ids), inter, len(contig_ids - kmer_ids)))
        print(f"k{k}: kmers={len(kmer_ids)} contigs={len(contig_ids)} "
              f"intersection={inter} kmers_seuls={len(kmer_ids - contig_ids)} "
              f"contigs_seuls={len(contig_ids - kmer_ids)}")

    # Tableau d'effectifs.
    csv_path = os.path.join(outdir, "venn_counts.csv")
    with open(csv_path, "w") as fh:
        fh.write("k,kmers_seuls,intersection,contigs_seuls\n")
        for k, ko, inter, co in rows:
            fh.write(f"{k},{ko},{inter},{co}\n")
    print(f"  -> {csv_path}")

    # Un diagramme par k.
    for k in ks:
        kmer_ids, contig_ids = ids[k]
        fig, ax = plt.subplots(figsize=(5, 5))
        venn_pair(ax, kmer_ids, contig_ids, k, show_title=False)
        fig.suptitle(f"Venn kmers vs contigs - entropie k{k}")
        for ext in ("png", "pdf"):
            out = os.path.join(outdir, f"venn_k{k}.{ext}")
            fig.savefig(out, dpi=200, bbox_inches="tight")
            print(f"  -> {out}")
        plt.close(fig)

    # Panneau recapitulatif (grille la plus carree possible).
    ncols = min(len(ks), 2)
    nrows = math.ceil(len(ks) / ncols)
    fig, axes = plt.subplots(nrows, ncols, figsize=(5 * ncols, 5 * nrows))
    axes = list(axes.flat) if hasattr(axes, "flat") else [axes]
    for ax, k in zip(axes, ks):
        venn_pair(ax, ids[k][0], ids[k][1], k)
    for ax in axes[len(ks):]:        # masque les cases vides eventuelles
        ax.axis("off")
    fig.suptitle("Venn kmers vs contigs par filtre d'entropie")
    fig.tight_layout()
    for ext in ("png", "pdf"):
        out = os.path.join(outdir, f"venn_all.{ext}")
        fig.savefig(out, dpi=200, bbox_inches="tight")
        print(f"  -> {out}")
    plt.close(fig)


def main(argv=None):
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--input", default=os.path.expanduser("~/complexity"),
                        help="dossier contenant les fichiers .fa (defaut: ~/complexity)")
    parser.add_argument("--outdir", default=None,
                        help="dossier de sortie (defaut: <input>/venn)")
    args = parser.parse_args(argv)

    input_dir = os.path.expanduser(args.input)
    if not os.path.isdir(input_dir):
        sys.exit(f"Erreur: dossier introuvable: {input_dir}")
    outdir = args.outdir or os.path.join(input_dir, "venn")
    os.makedirs(outdir, exist_ok=True)

    process(input_dir, outdir)
    print("Termine.")


if __name__ == "__main__":
    main()
