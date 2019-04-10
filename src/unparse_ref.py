import os
import psutil
import sys
import bisect
from bitstring import BitArray
from Walvelettree import WalveletTree
import unparse_read as ur
import time

wt = []
lf = []
ind_lf = []
pos = ""

converter = {'A': BitArray('0b00'),'C': BitArray('0b01'),'G': BitArray('0b10'),'T': BitArray('0b11')}

def _invconvert_(b):
    for letter, val in converter.items():
        if val == b:
            return letter

def _getILF_(i):
    """
    Retourne le ième caractère de la left_colum
    """
    pus = 0
    ind = 0
    while pus <= (i-1):
        pus += lf[ind]
        ind += 1
    return ind-1

def _getIC_(i, char):
    """
    Selectionne le ième CHAR(A,C,G ou T) dans la left-column
    """
    global pos

    intval = int(converter[char].bin,2)
    if (lf[intval] == 0):
        print('caractère non présent')
        sys.exit
    if (lf[intval] < i ):
        print('Ième caractère trop grand')
    return ind_lf[intval]+i

def _lttobw_(i):
    """
    Retourne l'indice dans bw de l'élement en face du ième dans la lf
    """
    return i if (i < pos) else (-1 if (i == pos) else i-1)


def _prettyprinter_():
    """
    Pretty printer pour afficher les 2 colonnes (à gauche : lettres ordonnées, à droite : bw)
    """
    print('0 \t$ \t'+_invconvert_(wt.getC(0)))
    for i in range(1, len(wt.data)+1):
        if i < (pos):
            print(str(i)+"\t"+str(_getILF_(i))+'\t'+_invconvert_(wt.getC(i)))
        elif i > (pos):
            print(str(i)+"\t"+str(_getILF_(i))+'\t'+_invconvert_(wt.getC(i-1)))
        else :
            print(str(i)+"\t"+str(_getILF_(i))+'\t$')


def foundKmer(kmer):
    """
    Trouve les différentes occurences du kmer dans la séquence et renvoie True ou False si il existe au moins une occurence
    """
    n = len(kmer)-1
    char = kmer[-1]
    intval = int(converter[char].bin,2)
    # première lettre : on trouve la 1ère et dernière occurence de la lettre dans la colonne de gauche (ici, first et last)
    first, last = _getIC_(0, char), _getIC_(lf[intval]-1, char)
    # notfound = True
    n -= 1
    while (n >= 0):
        char = kmer[n]

        _first = _lttobw_(first)
        _last = _lttobw_(last)

        first = (first+1) if(_first == -1) else _first
        last = (last-1) if(_last == -1) else _last

        if (first > last):
            # print("Pas de match possible")
            return False

        # rank du char dans la transformée
        rank1 = wt.rank(converter[char], first-1)
        rank2 = wt.rank(converter[char], last)

        if (rank1 >= rank2):
            # print("Pas de match possible")
            return False
        first, last = _getIC_(rank1, char), _getIC_(rank2-1, char)

        n -=1

    return True

def parseBWC(handle):
    """
    Parse un fichier contenant une transformée de Burrows Wheeler compressée et le transforme en walvelet tree
    """
    global lf
    global wt
    global pos
    global ind_lf

    counter = 0
    lf = [0,0,0,0]
    ind_lf = [0,0,0,0]

    rac = BitArray()
    n1 = BitArray()
    n2 = BitArray()

    # Lis la position de fin
    n_bytes_pos = handle.read(1)
    n_pos = int.from_bytes(n_bytes_pos,"big")
    bytes_pos = []
    for i in range(n_pos):
        bytes_pos+= handle.read(1)
    pos = int.from_bytes(bytes_pos,"big")

    last_bytes = [handle.read(1),handle.read(1)]

    while True:
        c = handle.read(1)
        if not c:
            c_bytes = BitArray(last_bytes[0])
            if (last_bytes[1] == b'\x00'): # si pas de bourrage
                n = 4
            else:
                n = int.from_bytes(last_bytes[1],"big")
            for i in range (0,n):
                char = c_bytes[2*i:2*i+2]
                rac.append([char[0]])
                if c_bytes[2*i]:
                    n2.append([char[1]])
                else :
                    n1.append([char[1]])
                lf[int(char.bin,2)] += 1
            break


        c_bytes = BitArray(last_bytes[0])

        for i in range (0,4):
            char = c_bytes[2*i:2*i+2]
            rac.append([char[0]])
            if c_bytes[2*i]:
                n2.append([char[1]])
            else :
                n1.append([char[1]])
            lf[int(char.bin,2)] += 1

        last_bytes = [last_bytes[1],c]

    for i in range(4):
        ind_lf[i] = sum(lf[0:i])+1

    wt = WalveletTree(rac, n1, n2)

def main2(bwfilename, fastqfilename, k):
    """
    Même chose que main mais avec une barre de progression
    """
    global wt
    global converter

    with open(bwfilename, 'rb') as bw_file:
        parseBWC(bw_file)
    with open("".join(bwfilename.split('.')[0].split('_')[0:2])+"_f.fastq", 'w') as kmer_file:
        for kmerset in ur.transfer_kmer(fastqfilename, k):
            k_pres = BitArray()
            for kmer in kmerset:
                k_pres.append([foundKmer(kmer)])
            header = "@"
            header += ur.commun_read_ref(k_pres, k)
            header +=  ur.coverage_read_ref(k_pres, k)
            kmer_file.write(header)
    # with open(bwfilename.split('.')[0]+"_kmer.txt", 'w') as kmer_file:
    #     for kmerset in transfer_kmer(fastqfilename, k):
            # counterRead += 1
            # counter = 0
            # ProgressBar = None
            # for kmer in kmerset:
            #     counter += 1
            #     label = "Loading read n°"+str(counterRead)+" (" + str(counter) + "/" + str(len(kmerset)) + ")"
            #     if ProgressBar == None:
            #         ProgressBar = progressbar.ProgressBar(counter, len(kmerset), label = label, usePercentage = False)
            #     else:
            #         ProgressBar.updateProgress(counter, label)
            #     bornes = foundKmer(kmer)
            #     kmer_file.write(kmer+" "+str(bornes[0])+" "+str(bornes[1])+'\n')

def toFileBitArray(oufile, k_pres):
    """
    Prend un bit array et le write dans un fichier
    """
    outfile.write(k_pres.bytes)

def main(bwfilename, fastqfilename, k):
    """
    Prend un fichier sérialisé et un jeu de données de référence
    """
    global wt
    global converter
    process = psutil.Process(os.getpid())
    print(process.memory_info().rss)

    with open(bwfilename, 'rb') as bw_file:
        start_time = time.time()
        parseBWC(bw_file)
        print("ProcessBWT--- %s seconds ---" % (time.time() - start_time))
    print(process.memory_info().rss)
    with open("".join(bwfilename.split('.')[0].split('_')[0:2])+"_f.fastq", 'w') as kmer_file:
        for kmerset in ur.transfer_kmer(fastqfilename, k):
            k_pres = BitArray()
            for kmer in kmerset:
                k_pres.append([foundKmer(kmer)])
                print("Process1KMER--- %s seconds ---" % (time.time() - start_time))
            header = "@"
            header += ur.commun_read_ref(k_pres, k)
            header +=  ur.coverage_read_ref(k_pres, k)
            kmer_file.write(header)
            print("Process1READ--- %s seconds ---" % (time.time() - start_time))
            print(process.memory_info().rss)


if __name__ == "__main__":
    """
    arg1 : bw ref file
    arg2 : reads file
    arg3 : kmer size
    """
    main(sys.argv[1], sys.argv[2], int(sys.argv[3]))
