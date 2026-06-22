# Comparaisons des filtres d'entropie (kmers & contigs)

Deux scripts, qui comparent tous deux les **identifiants uniques de fusion**
extraits des fichiers FASTA du dossier `complexity` :

- `upset_entropy.py` : **UpSet plots** comparant les k entre eux (k3 vs k5 vs
  k8 vs k10), un plot pour les kmers et un pour les contigs.
- `venn_entropy.py` : **diagrammes de Venn** comparant, pour chaque k, les
  kmers vs les contigs.

## upset_entropy.py

Compare, entre les filtres d'entropie **k3, k5, k8 et k10**, les identifiants
uniques de fusion presents dans :

- les kmers (`kmers-entropyk*.fa`),
- les contigs (`contigs-entropyk*.fa`).

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

## venn_entropy.py

Pour **chaque k** (k3, k5, k8, k10), trace un diagramme de Venn comparant les
identifiants de fusion presents dans les kmers vs dans les contigs.

```bash
# Par defaut : lit ~/complexity et ecrit dans ~/complexity/venn
python venn_entropy.py

# Ou en precisant les dossiers
python venn_entropy.py --input /chemin/vers/complexity --outdir /chemin/sortie
```

Sorties (dans `--outdir`) :

| Fichier | Contenu |
|---|---|
| `venn_k{3,5,8,10}.png` / `.pdf` | un Venn kmers vs contigs par k |
| `venn_all.png` / `.pdf` | panneau recapitulatif (tous les k) |
| `venn_counts.csv` | effectifs (kmers_seuls, intersection, contigs_seuls) par k |

## kmer_counts.py

Construit un tableau du **nombre de kmers** associes a chaque identifiant unique
de fusion, pour chaque filtre k3/k5/k8/k10.

```bash
# Par defaut : lit ~/complexity et ecrit kmer_counts.csv dans le meme dossier
python kmer_counts.py

# Ou en precisant les dossiers
python kmer_counts.py --input /chemin/vers/complexity --outdir /chemin/sortie
```

Sortie : `kmer_counts.csv` avec une ligne par fusion et les colonnes
`fusion_id, k3, k5, k8, k10` (nombre de kmers, 0 si la fusion est absente du
filtre).

## Note de compatibilite

`upsetplot 0.9.0` necessite `pandas < 3.0` (voir `requirements.txt`). Le script
n'utilise pas `show_counts` d'upsetplot — incompatible avec `numpy >= 2` — et
ajoute les effectifs lui-meme, il fonctionne donc avec numpy 2.x.
