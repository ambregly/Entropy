# UpSet plots des filtres d'entropie (kmers & contigs)

`upset_entropy.py` compare, entre les filtres d'entropie **k3, k5, k8 et k10**,
les **identifiants uniques de fusion** presents dans les fichiers FASTA du
dossier `complexity`, et produit deux UpSet plots :

- un pour les kmers (`kmers-entropyk*.fa`),
- un pour les contigs (`contigs-entropyk*.fa`).

## Principe

Chaque en-tete FASTA a la forme :

```
>ELOVL5_chr6_53213614_53213640_PTP4A1_chr6_64286340_64286365_2744.kmer1
>OAZ1_chr19_2269717_2269743_ZBTB7A_chr19_4055220_4055245_2564.contig_1
```

La cle de comparaison est l'**identifiant unique de fusion**, c.-a-d. l'en-tete
prive du suffixe `.kmerN` / `.contig_N`. Une fusion est comptee comme presente
pour un k donne des qu'au moins un kmer (resp. contig) la concerne dans le
fichier correspondant.

## Installation

```bash
pip install -r requirements.txt
```

## Utilisation

```bash
# Par defaut : lit ~/complexity et ecrit dans ~/complexity/upset
python upset_entropy.py

# Ou en precisant les dossiers
python upset_entropy.py --input /chemin/vers/complexity --outdir /chemin/sortie
```

## Sorties (dans `--outdir`)

| Fichier | Contenu |
|---|---|
| `upset_kmers.png` / `.pdf` | UpSet plot des fusions par k (kmers) |
| `upset_contigs.png` / `.pdf` | UpSet plot des fusions par k (contigs) |
| `membership_kmers.csv` | tableau presence/absence (fusion x k) |
| `membership_contigs.csv` | idem pour les contigs |

## Note de compatibilite

`upsetplot 0.9.0` necessite `pandas < 3.0` (voir `requirements.txt`). Le script
n'utilise pas `show_counts` d'upsetplot — incompatible avec `numpy >= 2` — et
ajoute les effectifs lui-meme, il fonctionne donc avec numpy 2.x.
