# SV2018
Célestine SAUVAGE

# Filtrage de k-mers sans alignement sur référence

Ce TP est capable pour une séquence donnée et pour un fichier contenant des reads, de filtrer les k-mers communs à la référence et de les écire dans un fichier Fastq. Il a été rédigé en python3.

## How to ?
### Utilisation
(On suppose ici qu'on est à la racine du projet)
* Lancer l'environnement virtuel pour charger les librairies nécéssaires au bon fonctionnement du programme :

```bash
source venv/bin/activate
```

* Pour sérialiser la séquence de référence, faire :

**/!\\ PAS DE SÉQUENCE TROP GROSSE /!\\**
*car le programme prend (BEAUCOUP) trop de mémoire pour calculer le suffix array*

```bash
python3 src/parse.py data/<nom_du_fichier_de_la_référence>.fa <OPTION : numéro_du_contig>
```
ou

```bash
python3 src/parse.py data/<nom_du_fichier_de_la_référence>.fa <OPTION : numéro_du_contig>
```
Le fichier contenant la référence doit être au format *FASTA*.

* Pour procéder au filtrage des reads :

```bash
python3 src/unparse_ref.py data/<binaire_généré> data/<nom_du_fichier_des_reads>.fastq <taille_des_kmers>
```
Le fichier de reads devra être au format *FASTQ*.

### Exemple
Pour le fichier Test.fa contenant le séquence **ACGTTGCAAACCCTTTGGGC**, je lance la sérialisation avec :

```bash
python3 src/parse.py data/Test.fa
```

On obtient un fichier binaire `data/Test_0_bw` qui contient la transformée de Burrows-Wheeler
(voir partie Méthode).

Ensuite : la commande

```bash
python3 src/unparse_ref.py data/Test_0_bw data/Test_reads.fastq 3
```
Renvoie le fichier `Test_reads_f.fastq` contenant les reads filtrés.

## Méthode

*Comment le programme fait-il pour calculer les k-mers communs ?*

### Partie parsing/sérialization
Dans un premier temps, le fichier `parse.py` va calculer les indices du [*suffix_array*](https://en.wikipedia.org/wiki/Suffix_array) de la séquence. Puis en déduire [*la transformée de Burrows-Wheeler*](https://en.wikipedia.org/wiki/Burrows%E2%80%93Wheeler_transform), qui va ensuite être bytée et enregistrée dans un fichier binaire.

Pour reprendre notre exemple, nous avons notre chaîne **ACGTTGCAAACCCTTTGGGC**.
Le suffix_array sera :
```
20  $
7   AAACCCTTTGGGC
8   AACCCTTTGGGC
9   ACCCTTTGGGC
0   ACGTTGCAAACCCTTTGGGC
19  C
6   CAAACCCTTTGGGC
10  CCCTTTGGGC
11  CCTTTGGGC
1   CGTTGCAAACCCTTTGGGC
12  CTTTGGGC
18  GC
5   GCAAACCCTTTGGGC
17  GGC
16  GGGC
2   GTTGCAAACCCTTTGGGC
4   TGCAAACCCTTTGGGC
15  TGGGC
3   TTGCAAACCCTTTGGGC
14  TTGGGC
13  TTTGGGC
```
La transformée BW associée sera `CCAA$GGACACGTGTCTTGTC`
Pour serialiser notre transformée, on va associer à chaque caractère 2 bits tel que:
```python3
A : 00
C : 01
G : 10
T : 11
```
et stocker la position du `$` dans un entier.

Le fichier binaire généré aura cette forme-là :
```
<nombre_de_bytes_pour_la_position_$> [bytes_position] [char1,char2,char3,char4] [char5,char6,char7,char8]... [char_n, bourrage] [nombre_de_bourrage]
```

Dans notre cas:
```bash
# lecture little endian
0000000 0401 a150 b71b 00ed
```

### Partie indexation
On va ensuite désérialiser la transformée et la stocker dans un [*wavelet tree*](https://en.wikipedia.org/wiki/Wavelet_Tree) en reprenant la même représentation pour les caractères:
```
racine : A,C -> 0
         G,T -> 1
gauche : A -> 0
         C -> 1
droite : G -> 0
         T -> 1
```
Ce qui nous donnera le wavelet tree suivant pour notre exemple
```
      00001100001111011110
1100010111            0001011101
```
On va en plus stocker le nombre d'occurences de chaque caractère dans une liste et la position du `$` pour pouvoir utiliser la méthode du [*FM_index*](https://en.wikipedia.org/wiki/FM-index).

### Partie kmérisation
Une fois notre structure mise en place pour la séquence de référence, on va lire les reads 1 à 1 et construire leur *set* de kmer associé. Pour chaque kmer, on va stocker aussi ses indices d'occurence pour pourvoir, a posteriori, calculer le `%` de coverage.

### Partie filtrage / coverage
Dans notre exemple: le filtrage/ calcul de coverage sera:
```
@kmers filtered  : ACG, CGT, GTT, TTG, TGC,
@5 / 10
@58.33
ACGTTGCGGAGA

@kmers filtered  : CGT,
@2 / 4
@100.00
CGTCGT

@kmers filtered  :
@0 / 4
@0.00
TTCCGG
```

## Implémentation, complexité en temps et en espace
Après avoir expliqué les méthodes utilisées pour l'indexation des kmers, nous allons rentrer plus en détail dans la partie implémentation. Nous allons supposer que les caractères sont réellement codés sur 2 bits.

### Libs utilisées
Pour faciliter l'implémentation du code, j'ai utilisé le module Python [**bitstring**](https://pythonhosted.org/bitstring/) qui permet la création et manipulation de structures binaire.

### Suffix_array
#### Implémentation
[*@source code*](https://www.geeksforgeeks.org/suffix-array-set-2-a-nlognlogn-algorithm/)
[*@source papier*](https://web.stanford.edu/class/cs97si/suffix-array.pdf)

Pour l'implémentation, j'ai trouvé une méthode qui permet de calculer un suffix array en utilisant le fait qu'on essaye de trier des strings qui sont des suffixes d'une seule chaîne.
On trie tout d'abord les suffixes par rapport à leur 2 premiers caractères en les remplaçant par un rang (ici on compare au code ASCII de la lettre 'A').
On va donc trier une première fois les suffixes par rapport à ces 2 rangs.
Après ça, on assigne un nouveau rang à tous les suffixes, en fonction de leur 2 rangs précédents et l'autre rang en fonction du rang du suffixe i +2.
Etc,etc...

**/!\\** Je n'ai pas réussi à me libérer de la contrainte <Objet suffixe> alors que j'aurai pu coder cette algorithm seulement avec 4 listes de taille n d'entiers. Il y a également une liste en plus que je pourrais supprimer (la liste suffixe_array de fin).
Une des solutions aurait pu être de calculer des suffixes array intermédiaire pour chaque contigs.
#### Complexité
**Temps** O(n+ logn*(2n+nlogn))

**Espace** 4 listes d'entier de taille n.

### BWT
#### Implémentation
En Python, pas de difficulté : pour chaque indice du tableau, faire `T[si-1]`
#### Complexité
**Temps** O(n)

**Espace** Dans l'idéal : `2b*(n-1)`

### Taille du fichier sérialisé
#### Implémentation
On met la BWT dans un fichier binaire en écrivant un byte pour 4caractères ensemble.
**/!\\** une première méthode que j'ai faite consistait à remettre à les caractères (pas en binaire) et en remplaçant les répétitions par des `<nb_caractères><caractère>`, qui était plus rapide à sérialiser, désérialiser, mais le fichier prenait plus de place.

#### Complexité
**Temps** O(n)

**Espace** Dans l'idéal : `2b*(n-1) + 8b + min(8,log2(position_$)) + bourrage`

### Wavelet tree
#### Implémentation
Comme les caractères sont déjà sous forme binaire dans le fichier sérialisé, il suffit de les lire 1 à 1 (en stockant les 2 derniers lus pour la gestion des caractères de bourrage) et de stocker le bit[0] dans le bit array root et de stocker dans le sous-arbre gauche ou droite le bit[1].

Pour ce qui est du *rank* des caractères, on le fait en temps constant en shiftant les bits et en calculant le [*Hamming weight*](https://en.wikipedia.org/wiki/Hamming_weight)(nombre de bits à 1) dans le root et le bon sous-arbre.
#### Complexité
**Temps** Construction -> O(n), rank -> O(1)

**Espace** Dans l'idéal : `2b*(n-1)`

### FM method
#### Implémentation
En plus du walvelet tree, on va stocker le nombre d'occurence de chaque nucléotide (liste de 4 entiers), où apparaît la première occurence du mot + la position du `$` dans la BWT.
Pour calculer la position du ième C dans la colonne de gauche, on va récupérer la 1ère occurence de la lettre et ajouter i.
#### Complexité
**Temps** O(1)

**Espace** walvelet_tree + 9 entiers.

### Kmérisation
#### Implémentation
On lit les caractères un à un pour construire un kmer, si il a déjà été stocké, on ajoute sa position occurence dans la liste, sinon on l'ajoute dans la liste des kmers + position d'occurence.
#### Complexité
**Temps** O(n)

**Espace** `(n-k)+1` kmers dans le pire des cas + 2n entiers  

### Calcul kmers communs / coverage
#### Implémentation
Après avoir filtré le read sur la séquence, on obtient un vecteur de bits pour savoir si le i_ème kmer et oui ou non dans la séquence.

**Kmers communs** : On compte pour chaque kmer présent son nombre d'occurences.
**Coverage** : On compte pour chaque occurence de chaque kmer, les ensembles du read qu'il couvre et on calcule l'union de ses ensembles pour avoir la couverture totale.

#### Complexité
**Temps kmers communs** O(n)
**Temps coverage** `O(n * 2log2(n))` (pour construire les ensembles trié par ordre croissant) + `O(2n)` pour l'union des ensembles.

## Benchmarks
### Calcul pour le 1er contig du suffix array + sérialisation
```
ParseFile--- 0.011195659637451172 seconds ---
35676160 octets
Traitement de la référence de longueur 2021796
Suffix_array--- 106.47510600090027 seconds ---
445566976 octets
BWT--- 107.09699606895447 seconds ---
445566976 octets
	 BW compressée
BWT compressée--- 127.55657362937927 seconds ---
445566976 octets
```
### Calcul pour le 10ème contig du suffix array + sérialisation
```
python3 src/parse.py data/Acanthamoeba.fa 10
ParseFile--- 0.07288146018981934 seconds ---
35299328 octets
Traitement de la référence de longueur 722550
Suffix_array--- 34.089083433151245 seconds ---
179957760 octets
BWT--- 34.26090598106384 seconds ---
179957760 octets
	 BW compressée
BWT compressée--- 41.64849781990051 seconds ---
179957760 octets
```

### Désérialisation + walvelet tree pour le 1er contig
```
13455360
ProcessBWT--- 116.37714982032776 seconds ---
14249984
```

### Indexation des kmers sur le 1er contig avec k = 4, taille read = 45
```
Process1KMER--- 117.04473853111267 seconds ---
Process1KMER--- 117.57030701637268 seconds ---
Process1KMER--- 118.03926587104797 seconds ---
Process1KMER--- 118.50457286834717 seconds ---
Process1KMER--- 119.05024075508118 seconds ---
Process1KMER--- 119.64877104759216 seconds ---
Process1KMER--- 120.29394745826721 seconds ---
Process1KMER--- 120.94233322143555 seconds ---
Process1KMER--- 121.53990650177002 seconds ---
Process1KMER--- 122.06502747535706 seconds ---
Process1KMER--- 122.61584544181824 seconds ---
Process1KMER--- 123.16900682449341 seconds ---
Process1KMER--- 123.69396018981934 seconds ---
Process1KMER--- 124.18233132362366 seconds ---
Process1KMER--- 124.6598973274231 seconds ---
Process1KMER--- 125.22231888771057 seconds ---
Process1KMER--- 125.69437384605408 seconds ---
Process1KMER--- 126.25497674942017 seconds ---
Process1KMER--- 126.96148562431335 seconds ---
Process1KMER--- 127.48319673538208 seconds ---
Process1KMER--- 128.0254967212677 seconds ---
Process1KMER--- 128.6562762260437 seconds ---
Process1KMER--- 129.1669201850891 seconds ---
Process1KMER--- 129.73000025749207 seconds ---
Process1KMER--- 130.22063517570496 seconds ---
Process1KMER--- 130.83118653297424 seconds ---
Process1KMER--- 131.42218112945557 seconds ---
Process1KMER--- 132.03872990608215 seconds ---
Process1KMER--- 132.60918807983398 seconds ---
Process1KMER--- 133.2208743095398 seconds ---
Process1KMER--- 133.8462929725647 seconds ---
Process1READ--- 133.84648180007935 seconds ---
17612800
```
~ - de 1 seconde par kmer + calcul instantané des kmers communs + taux de coverage

### Indexation des kmers sur le 1er contig avec k = 20, taille read = 45

```
Process1KMER--- 121.8738603591919 seconds ---
Process1KMER--- 124.11825442314148 seconds ---
Process1KMER--- 126.44388055801392 seconds ---
Process1KMER--- 129.0914342403412 seconds ---
Process1KMER--- 131.3854010105133 seconds ---
Process1KMER--- 133.73787879943848 seconds ---
Process1KMER--- 135.88102650642395 seconds ---
Process1KMER--- 138.27094841003418 seconds ---
Process1KMER--- 140.45185256004333 seconds ---
Process1KMER--- 142.84943652153015 seconds ---
Process1KMER--- 145.4413845539093 seconds ---
Process1KMER--- 147.40780425071716 seconds ---
Process1KMER--- 149.4126353263855 seconds ---
Process1KMER--- 151.24889373779297 seconds ---
Process1KMER--- 153.17265224456787 seconds ---
Process1KMER--- 155.27994894981384 seconds ---
Process1KMER--- 157.65046954154968 seconds ---
Process1KMER--- 159.80158114433289 seconds ---
Process1KMER--- 161.9048662185669 seconds ---
Process1KMER--- 163.83560490608215 seconds ---
Process1KMER--- 165.74338698387146 seconds ---
Process1KMER--- 167.91950225830078 seconds ---
Process1KMER--- 170.25797176361084 seconds ---
Process1KMER--- 172.46156573295593 seconds ---
Process1READ--- 172.46169233322144 seconds ---
17616896
```
