from bitstring import BitArray

class WalveletTree:

    def __init__(self, data, l, r):
        self.left = l
        self.right = r
        self.data = data
        # self.converter = {'A': BitArray('0b00'),'C': BitArray('0b01'),'G': BitArray('0b10'),'T': BitArray('0b11')}
    #
    def getC(self, i):
        """
        Retourne le caractère situé à l'indice i
        """
        b =self.data[i]
        f_val = self._rankbitdata_(b, i, self.data)
        if b:
            b2 = self.right[f_val]
        else:
            b2 = self.left[f_val]
        return BitArray(list([b,b2]))


    def rank(self, c_val, i):
        """
        Renvoie le nombre d'occurrence du caractère c_val (sous forme binaire), de 0 à l'indice i
        en utilisant le décalage de bits et la méthode count
        """
        i += 1
        assert(i <= len(self.data))
        if c_val[0] : # G/T
            val_f = (self.data >> (len(self.data) - i)).count(True)
            if c_val[1] : #G
                val_s = (self.right >> (len(self.right) - val_f)).count(True)
            else : #T
                val_s = val_f - (self.right >> (len(self.right) - val_f)).count(True)
        else : # A/C
            val_f = i - (self.data >> (len(self.data) - i)).count(True)
            if c_val[1] : #C
                val_s = (self.left >> (len(self.left) - val_f)).count(True)
            else : #A
                val_s = val_f - (self.left >> (len(self.left) - val_f)).count(True)
        return val_s

    def _rankbitdata_(self, b, i, data):
        if b:
            return (data >> (len(data) - i)).count(True)
        else :
            return i - (data >> (len(data) - i)).count(True)

    def select(self, c_val, i):
        """
        Renvoie la ième lettre c_val dans le chaîne
        """
        if c_val[0]:
            val_s = self.select_dicho(c_val[1], i, self.right)
        else :
            val_s = self.select_dicho(c_val[1], i, self.left)

        return self.select_dicho(c_val[0], val_s, self.data) -1


    def select_dicho(self, b, nbO, data):
        l, r = nbO, len(data)
        # assert(self._rankbitdata_(b, r, data) > i)
        # val_rank = self._rankbitdata_(b, nbO, data)
        while l < r:
            k = (l+r)//2
            val_rank = self._rankbitdata_(b, k, data)
            if (val_rank < nbO):
                l = k + 1
            else :
                r = k
        return l
