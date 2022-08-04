# ####_______________ Parsing Cigar String : from SAM format
# M=0 #M  BAM_CMATCH      0
# I=1 #I  BAM_CINS        1
# D=2 #D  BAM_CDEL        2
# N=3 #N  BAM_CREF_SKIP   3
# S=4 #S  BAM_CSOFT_CLIP  4
# H=5 #H  BAM_CHARD_CLIP  5
# P=6 #P  BAM_CPAD        6
# E=7 #=  BAM_CEQUAL      7
# X=8 #X  BAM_CDIFF       8
# B=9 #B  BAM_CBACK       9
# NM=10 #NM       NM tag  10
# conRef  =       [M, D, N, E, X] # these ones "consume" the reference
# conQuery=       [M, I, S, E, X] # these ones "consume" the query
# conAln  =       [M, I, D, N, S, E, X] # these ones "consume" the alignments
# ####________________
cigar_dict = { 0 : 'M' , 1 : 'I' , 2 : 'D' , 3 : 'N' , 4 : 'S' , 5 : 'H' , 6 : 'P' , 7 : 'E', 8 : 'X' , 9 : 'B' , 10 : 'NM'}

class CigarOperation:
    cigar_dict = { 0 : 'M' , 1 : 'I' , 2 : 'D' , 3 : 'N' , 4 : 'S' , 5 : 'H' , 6 : 'P' , 7 : 'E', 8 : 'X' , 9 : 'B' , 10 : 'NM'}
    def __init__(self, match = 1, mismatch = 9, gap_open = [16, 41] , gap_extend = [2,1] ):
        self.match = match
        self.mismatch = mismatch
        self.gap_open = gap_open 
        self.gap_extend = gap_extend
    def operation(self, opt, l):
        '''get a given operation from cigar string.'''
        self.l = int(l)
        opt = 'M' if opt == '=' else opt #turn operation into a M if it's a =
        opt = opt if isinstance(opt, str) else cigar_dict[opt]
        return getattr(self, 'opt_' + str(opt), "Non supported character")()
    def opt_M(self):
        '''get score from match'''
        return self.match * self.l 
    def opt_D(self):
        '''score of deletion'''
        return -1 * min(self.gap_open[0] + self.l * self.gap_extend[0]  , self.gap_open[1] + self.l * self.gap_extend[1])
    def opt_I(self):
        '''insertion : same operation as deletion.'''
        return CigarOperation.opt_D(self)
    def opt_N(self):
        '''skip sequence'''
        return 0
    def opt_S(self):
        '''soft_clipping: Skip'''
        return 0
    def opt_H(self):
        '''hard_clipping: Skip'''
        return 0
    def opt_P(self):
        '''Padding : Skip'''
        return 0
    def opt_E(self):
        '''equal'''
        return CigarOperation.opt_M(self)
    def opt_X(self):
        '''mismatch'''
        return -1 * self.mismatch * self.l


