import sys
import gzip
import bisect

kmer_set = []
kmer_ind = []
size_read = 0

def _sum(ens):
    size = 0
    for inter in ens:
        size += inter[1]-inter[0]
    return size

def _union(l_d, l_f):
    """
    Union d'ensembles
    """
    ens_f = []
    nb_s = 0
    d_courant = 0
    while l_d:
        if l_d[0] < l_f[0]:
            pos_debut = l_d.pop(0)
            if nb_s == 0:
                d_courant = pos_debut
            nb_s += 1
        elif l_d[0] > l_f[0]:

            pos_fin = l_f.pop(0)
            nb_s -= 1
            if nb_s == 0:
                n_int = (d_courant, pos_fin)
                ens_f.append(n_int)
        else:
            l_d.pop(0)
            l_f.pop(0)

    if l_f:
        pos_fin = l_f[-1]
        n_int = (d_courant, pos_fin)
        ens_f.append(n_int)

    return ens_f


def coverage_read_ref(bit_vector, k):
    """
    Renvoie le pourcentage de couverture des kmers en fonction de la taille du read
    """
    global size_read
    global read

    l_d = [] #debut des ensembles
    l_f = [] #fin des ensembles
    for i, b in enumerate(bit_vector):
        if b:
            for occ in kmer_ind[i]:
                bisect.insort(l_d, occ-k)
                bisect.insort(l_f, occ)
    union = _union(l_d,l_f)
    return "\n@"+ ("%.2f" % (_sum(union)*100 / size_read)) +"\n"+read+"\n"

def commun_read_ref(bit_vector, k):
    """
    Renvoie le nombre de kmer communs
    """
    global size_read
    global read

    nbkmer_tot = size_read - k+1
    nbkmer_com = 0
    kmers = ''
    for i, b in enumerate(bit_vector):
        if b:
            focc = kmer_ind[i][0]
            kmers += read[focc-k : focc]+", "
            occ = len(kmer_ind[i])
            nbkmer_com += occ
    return "kmers filtered  : "+kmers +"\n@"+str(nbkmer_com) +" / "+ str(nbkmer_tot)

def build_kmer(k):
    """
    Construit un ensemble unique des kmers du read
    Stocke les indexs de chaque kmer
    """
    global kmer_set #stocke l'ensemble des kmers
    global kmer_ind #stocker les indices associés
    global size_read
    global read

    kmer_set = dict()
    kmer_ind = []
    ind = 0

    kmer = read[0:k] # premier kmer
    yield kmer
    size_read = len(read) # taille du read

    kmer_set[kmer] = 0 # je stocke l'indice du premier kmer dans le tableau
    kmer_ind.append([k]) #je stocke l'indice de fin de kmer
    ind = 1 #indice pour le tableau
    while k < size_read:
        kmer = kmer[1:]+read[k] # je lis les kmers
        k += 1
        if kmer in kmer_set: # si il y est déjà
            kmer_ind[kmer_set[kmer]] += [k] #
        else :
            kmer_set[kmer] = ind
            kmer_ind.append([k])
            ind += 1
            yield kmer


def fastq_parser(handle):
    """
    Parse un fichier fastq
    """
    lines = []

    for line in handle:
        if line[0] == '@':
            lines = []
            continue
        if line[0] == '+':
            yield ''.join(lines).replace(" ", "").replace("\r", "").replace("N","")
        else:
            lines.append(line.rstrip())

    yield ''.join(lines).replace(" ", "").replace("\r", "").replace("N","")


def transfer_kmer(fastqfilename, k):
    """
    Transfert les kmers
    """
    global read
    with open(fastqfilename, 'r') as fastq_file:
        parsed_fasta = fastq_parser(fastq_file)
        for read in parsed_fasta:
            yield build_kmer(k)

if __name__ == "__main__":
    # execute only if run as a script
    main(sys.argv[1])
