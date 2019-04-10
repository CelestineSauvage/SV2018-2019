import time
import sys
import bisect
import os
import psutil
from bitstring import BitArray
from math import log
import numpy as np

suffix_array = []
bw = []
compbw = []
seq = []
pos = -1

converter = {'A': BitArray('0b00'),'C': BitArray('0b01'),'G': BitArray('0b10'),'T': BitArray('0b11')}


def bytes_needed(n):
    if n == 0:
        return 1
    return int(log(n, 256)) + 1

class mySuffix:

    __slots__ = ['index', 'rank']

    def __lt__(a,b):
        return (True if a.rank[1] < b.rank[1] else False) if (a.rank[0] == b.rank[0]) else (True if a.rank[0] < b.rank[0] else False)

    def __gt__(a,b):
        return (True if a.rank[1] >= b.rank[1] else False) if (a.rank[0] == b.rank[0]) else (True if a.rank[0] >= b.rank[0] else False)

    def __init__(self, index):
        self.index = index
        self.rank = [0,0]

    def toString(self):
        return (str(self.index)+'\t'+str(self.rank[0])+'\t'+str(self.rank[1]))

def toPrint(tab):
    for i in range(len(tab)):
        print(tab[i].toString())
    print("------------------------")

def suffix_array_rank():
    global suffix_array

    n = len(seq)
    suffixes= np.array([None] * n)

    for i in range(n-1):
        suffixes[i] = mySuffix(i)
        suffixes[i].rank[0] = ord(seq[i]) - ord('A')
        suffixes[i].rank[1] = ord(seq[i+1]) - ord('A') #if ((i+1) < n) else -1
    suffixes[n-1] = mySuffix(n-1)
    suffixes[n-1].rank[0] = ord(seq[n-1]) - ord('A')
    suffixes[n-1].rank[1] = -1

    suffixes = sorted(suffixes)
    ind = np.array([0] * n)

    k = 4
    while (k < 2*n):
        label = "Suffix array (" + str(k) + "/" + str(n) + ")"

        rank = 0
        prev_rank = suffixes[0].rank[0]
        suffixes[0].rank[0] = rank;
        ind[suffixes[0].index] = 0;

        for i in range(n):
            if ((suffixes[i].rank[0] == prev_rank) and
                    (suffixes[i].rank[1] == suffixes[i-1].rank[1])):
                prev_rank = suffixes[i].rank[0]
                suffixes[i].rank[0] = rank
            else :
                prev_rank = suffixes[i].rank[0]
                rank += 1
                suffixes[i].rank[0] = rank
            ind[suffixes[i].index] = i

        for i in range(n):
            next_index = suffixes[i].index + k//2
            suffixes[i].rank[1] =   suffixes[ind[next_index]].rank[0] if (next_index < n)  else -1

        suffixes = sorted(suffixes)

        k *= 2

    suffix_array = [n]
    for i in range(n):
        suffix_array.append(suffixes[i].index)

def BWT():
    global bw
    global seq
    global pos

    mod = 1
    count = 0

    for si in suffix_array:
        if si == 0:
            pos = count
        else:
            bw.append(converter[seq[si-1]])
            count += 1

def compressBW(bw_file, size):
    """
    Compresse la BW dans un fichier binaire
    """
    global bwinv
    global bw
    global pos

    count = 0
    bourre = 0

    n_bytes_pos = bytes_needed(pos)
    bytes_pos = (pos).to_bytes(n_bytes_pos, byteorder="big")
    bw_file.write(bytes([n_bytes_pos]))
    bw_file.write(bytes_pos) # on écrit la position du $


    while count < size:
        try: # on prend 4 caractères par 4 et on l'écrit sous forme de byte
            to_bytes = BitArray()
            l = bw[count:count+4]
            to_bytes = to_bytes.join(l)
            bw_file.write(to_bytes.bytes)
            count += 4
        except: # on bourre avec des 00 les derniers caractères
            to_bytes = BitArray()
            l = bw[count::]
            n = len(l)
            bourre = len(l)
            while n < 4:
                l.append(BitArray('0b00'))
                n += 1
            to_bytes = to_bytes.join(l)
            bw_file.write(to_bytes.bytes)
            count += 4
    bw_file.write(bytes([bourre]))


def _converte_(line):
    return list(map(lambda x: converter[x], line))

def fastaParser(handle):
    """
    Parse fasta file
    """
    lines = []

    for line in handle:
        if line[0] == '>':
            break
        lines.append(line.rstrip())
        break

    for line in handle:
        if line[0] == '>':
            yield ''.join(lines).replace(" ", "").replace("\r", "").replace("N","")
            lines = []
            continue
        lines.append(line.rstrip())

    yield ''.join(lines).replace(" ", "").replace("\r", "").replace("N","")

def _printsa_(s):
    for val in suffix_array:
        print((val,s[val::]))

def main(fastafilename, num = -1):
    global suffix_array
    global bw
    global compbw
    global seq

    process = psutil.Process(os.getpid())
    print(process.memory_info().rss)

    with open(fastafilename, 'r') as fasta_file:
        start_time = time.time()

        parsed_fasta = fastaParser(fasta_file)

        seq = []
        counter = 0
        if num == -1:
            for part in parsed_fasta:
                seq += part

        for part in parsed_fasta:
            if counter == num :
                seq = part
                break
            counter += 1
        print("ParseFile--- %s seconds ---" % (time.time() - start_time))
        print(process.memory_info().rss)

    suffix_array = []
    bw = []
    compbw = []
    print('Traitement de la référence de longueur '+str(len(seq)))
    suffix_array_rank()
    print("Suffix_array--- %s seconds ---" % (time.time() - start_time))
    print(process.memory_info().rss)

    BWT()
    print("BWT--- %s seconds ---" % (time.time() - start_time))
    print(process.memory_info().rss)

    print("\t BW compressée")
    with open(fastafilename.split('.')[0]+"_"+str(counter)+"_bw", 'wb') as bw_file:
        compressBW(bw_file,len(bw))
            #     for val in compbw:
            #         if val[0] == 1:
            #             bw_file.write(val[1])
            #         else :
            #             bw_file.write(str(val[0])+ val[1])
            # counter += 1
    print("BWT compressée--- %s seconds ---" % (time.time() - start_time))
    print(process.memory_info().rss)


if __name__ == "__main__":
    # execute only if run as a script
    if len(sys.argv) == 3:
        main(sys.argv[1], num = int(sys.argv[2]))
    else :
        main(sys.argv[1])
