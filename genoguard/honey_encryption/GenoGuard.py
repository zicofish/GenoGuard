'''
Created on Jun 17, 2014

@author: zhihuang

This is the major implementation for the paper:
Z. Huang, E. Ayday, J. Fellay, J.-P. Hubaux and A. Juels. GenoGuard: Protecting Genomic Data against Brute-Force Attacks. 36th IEEE Symposium on Security and Privacy (S&P 2015), San Jose, CA, USA, 2015.
'''
import math
import random
import numpy as np
import gmpy2
from gmpy2 import mpz
import struct
from Crypto.Hash import HMAC
from Crypto.Cipher import AES
from Crypto import Random
import time


class BaseGenoGuard(object):
    PREC = 10**-4
    powOfThree = []
    def __init__(self, SNPRefList, maxSeed):
        assert (maxSeed + 1) >= mpz(3)**mpz(len(SNPRefList))
        self.SNPRefList = SNPRefList
        self.maxSeed = maxSeed
        
        while len(self.powOfThree) < len(SNPRefList):    # For efficiency, pre-compute large powers of 3, as they will be used heavily in encoding and decoding.
            self.powOfThree.append(mpz(3) ** mpz(len(self.powOfThree)))
            
    def getMaxSeed(self):
        return self.maxSeed
    
    def encode(self, seq):
        assert len(seq) == len(self.SNPRefList)
        
        L = [-1, -1]; U = [-1, -1]
        L[0] = mpz(0); U[0] = mpz(self.maxSeed)
        n = len(seq)
        for i in range(1, n + 1):
            probs = self.condProb(seq[:i])
            sortInd = np.argsort(probs)
            assignedRangeSize = [0] * 3
            available = U[(i-1)%2] - L[(i-1)%2] + 1
            if probs[sortInd[0]] * available < self.powOfThree[n-i]:
                assignedRangeSize[sortInd[0]] = self.powOfThree[n-i]
            else:
                assignedRangeSize[sortInd[0]] = mpz(gmpy2.ceil(probs[sortInd[0]] * available))
            available -= assignedRangeSize[sortInd[0]]
            probs[sortInd[0]] = 0
            probs = probs / sum(probs)  #re-normalize the remaining two probabilities
            if probs[sortInd[1]] * available < self.powOfThree[n-i]:
                assignedRangeSize[sortInd[1]] = self.powOfThree[n-i]
            else:
                assignedRangeSize[sortInd[1]] = mpz(gmpy2.ceil(probs[sortInd[1]] * available))
            assignedRangeSize[sortInd[2]] = available - assignedRangeSize[sortInd[1]]
            L[i%2] = L[(i-1)%2] + sum(assignedRangeSize[:seq[i-1]])
            U[i%2] = L[(i-1)%2] + sum(assignedRangeSize[:(seq[i-1]+1)]) - 1
        #Don't use random.randint for big integer. Sometimes they might be correct on a system, but sometimes they are not. 
        #I encounter this problem when running this script with condor on icsil1-cluster.epfl.ch. The function generates very strange and unexpected integers.
        #seed = random.randint(L[n%2], U[n%2])
        tmp = random.randint(0, 2**30)
        randomState = gmpy2.random_state(tmp)
        seed = gmpy2.mpz_random(randomState, U[n%2] - L[n%2] + 1) + L[n%2]
        return seed
    
    def decode(self, seed):
        L = [-1, -1]; U = [-1, -1]
        L[0] = mpz(0); U[0] = mpz(self.maxSeed)
        n = len(self.SNPRefList)
        seq = [-1] * n
        for i in range(1, n + 1):
            probs = self.condProb(seq[:i])
            sortInd = np.argsort(probs)
            assignedRangeSize = [0] * 3
            available = U[(i-1)%2] - L[(i-1)%2] + 1
            if probs[sortInd[0]] * available < self.powOfThree[n-i]:
                assignedRangeSize[sortInd[0]] = self.powOfThree[n-i]
            else:
                assignedRangeSize[sortInd[0]] = mpz(gmpy2.ceil(probs[sortInd[0]] * available))
            available -= assignedRangeSize[sortInd[0]]
            probs[sortInd[0]] = 0
            probs = probs / sum(probs)  #re-normalize the remaining two probabilities
            if probs[sortInd[1]] * available < self.powOfThree[n-i]:
                assignedRangeSize[sortInd[1]] = self.powOfThree[n-i]
            else:
                assignedRangeSize[sortInd[1]] = mpz(gmpy2.ceil(probs[sortInd[1]] * available))
            assignedRangeSize[sortInd[2]] = available - assignedRangeSize[sortInd[1]]
            if seed < L[(i-1)%2] + sum(assignedRangeSize[:1]):
                seq[i-1] = 0
            elif seed < L[(i-1)%2] + sum(assignedRangeSize[:2]):
                seq[i-1] = 1
            elif seed < L[(i-1)%2] + sum(assignedRangeSize[:3]):
                seq[i-1] = 2
            else:
                raise ValueError("Invalid Seed")
            L[i%2] = L[(i-1)%2] + sum(assignedRangeSize[:seq[i-1]])
            U[i%2] = L[(i-1)%2] + sum(assignedRangeSize[:(seq[i-1] + 1)]) - 1
        return seq  
    
    def condProb(self, prefixSeq):
        return 0
    
    '''
    Pack a big integer into binary dat.
    '''
    def int_to_packed(self, val, bitLen):
        assert bitLen % 32 == 0
        num_words = bitLen / 32
        max_int = 2 ** mpz(bitLen) - 1
        max_word = 2 ** mpz(32) - 1
        words = []
        for _ in range(num_words):
            word = val & max_word
            words.append(word)
            val >>= 32
        words.reverse()
        return struct.pack('>%d%s' % (num_words, 'I'), *words)
    
    '''
    Unpack binary dat to a big integer.
    '''
    def packed_to_int(self, packed, bitLen):
        assert bitLen % 32 == 0
        num_words = bitLen / 32
        words = list(struct.unpack('>%d%s' % (num_words, 'I'), packed))
        words.reverse()
        val = mpz(0)
        for i, num in enumerate(words):
            word = mpz(num)
            word = word << (32 * i)
            val = val | word
        return val
    
    '''
    Password-based encryption.
    '''
    def pbeEncrypt(self, msg, password):
        assert len(msg) % 16 == 0
        S = Random.new().read(8) #64-bit random salt
        h = HMAC.new(password)
        h.update(S)
        derived_key = h.hexdigest()[:16]
        
        iv = Random.new().read(AES.block_size)
        cipher = AES.new(derived_key, AES.MODE_CBC, iv)
        ciphertext = cipher.encrypt(msg)
        assert len(ciphertext) == len(msg)
        return (S, iv, ciphertext)
    
    '''
    Password-based decryption.
    '''
    def pbeDecrypt(self, ext_c, password):
        S, iv, ciphertext = ext_c
        h = HMAC.new(password)
        h.update(S)
        derived_key = h.hexdigest()[:16]
        
        decipher = AES.new(derived_key, AES.MODE_CBC, iv)
        msg = decipher.decrypt(ciphertext)
        assert len(ciphertext) == len(msg)
        return msg
    
    '''
    Honey encryption including encoding and password-based encryption.
    '''
    def hgEncrypt(self, seq, password):
        seed = self.encode(seq)
        msg = self.int_to_packed(seed, bitLen)
        ext_c = self.pbeEncrypt(msg, password)
        return ext_c
    
    '''
    Honey decryption including password-based decryption and decoding.
    '''
    def hgDecrypt(self, ext_c, password):
        msg = self.pbeDecrypt(ext_c, password)
        seed = self.packed_to_int(msg, bitLen)
        seq = self.decode(seed)
        return seq
    
    '''
    Just some stupid code to calculate the running time of each function (for producing Figure 13 in the paper)
    '''
    def perfTest(self, seq, password):
        encode_start = time.time()
        seed = self.encode(seq)
        encode_end = time.time()
        msg = self.int_to_packed(seed, self.bitLen)
        pbe_enc_start = time.time()
        ext_c = self.pbeEncrypt(msg, password)
        pbe_enc_end = time.time()
        
        pbe_dec_start = time.time()
        msg = self.pbeDecrypt(ext_c, password)
        pbe_dec_end = time.time()
        seed = self.packed_to_int(msg, self.bitLen)
        decode_start = time.time()
        seq = self.decode(seed)
        decode_end = time.time()
        
        return (encode_end - encode_start, decode_end - decode_start, pbe_enc_end - pbe_enc_start, pbe_dec_end - pbe_dec_start)
    
    '''
    Compute the size of the interval that a sequence is encoded into
    Almost the same as encode(sefl, seq), but this time, we want the size of the final interval, but not a random integer in it.
    '''
    def seqIntervalSize(self, seq):
        assert len(seq) == len(self.SNPRefList)
        
        L = [-1, -1]; U = [-1, -1]
        L[0] = mpz(0); U[0] = mpz(self.maxSeed)
        n = len(seq)
        for i in range(1, n + 1):
            probs = self.condProb(seq[:i])
            sortInd = np.argsort(probs)
            assignedRangeSize = [0] * 3
            available = U[(i-1)%2] - L[(i-1)%2] + 1
            if probs[sortInd[0]] * available < self.powOfThree[n-i]:
                assignedRangeSize[sortInd[0]] = self.powOfThree[n-i]
            else:
                assignedRangeSize[sortInd[0]] = mpz(gmpy2.ceil(probs[sortInd[0]] * available))
            available -= assignedRangeSize[sortInd[0]]
            probs[sortInd[0]] = 0
            probs = probs / sum(probs)  #re-normalize the remaining two probabilities
            if probs[sortInd[1]] * available < self.powOfThree[n-i]:
                assignedRangeSize[sortInd[1]] = self.powOfThree[n-i]
            else:
                assignedRangeSize[sortInd[1]] = mpz(gmpy2.ceil(probs[sortInd[1]] * available))
            assignedRangeSize[sortInd[2]] = available - assignedRangeSize[sortInd[1]]
            L[i%2] = L[(i-1)%2] + sum(assignedRangeSize[:seq[i-1]])
            U[i%2] = L[(i-1)%2] + sum(assignedRangeSize[:(seq[i-1]+1)]) - 1
        return U[n%2] - L[n%2] + 1
    
    '''
    Compute the size of the interval that a seed falls into
    Almost the same as decode(self, seed), but this time, we want the size of the final interval, but not the original sequence.
    '''
    def seedIntervalSize(self, seed):
        L = [-1, -1]; U = [-1, -1]
        L[0] = mpz(0); U[0] = mpz(self.maxSeed)
        n = len(self.SNPRefList)
        seq = [-1] * n
        for i in range(1, n + 1):
            probs = self.condProb(seq[:i])
            sortInd = np.argsort(probs)
            assignedRangeSize = [0] * 3
            available = U[(i-1)%2] - L[(i-1)%2] + 1
            if probs[sortInd[0]] * available < self.powOfThree[n-i]:
                assignedRangeSize[sortInd[0]] = self.powOfThree[n-i]
            else:
                assignedRangeSize[sortInd[0]] = mpz(gmpy2.ceil(probs[sortInd[0]] * available))
            available -= assignedRangeSize[sortInd[0]]
            probs[sortInd[0]] = 0
            probs = probs / sum(probs)  #re-normalize the remaining two probabilities
            if probs[sortInd[1]] * available < self.powOfThree[n-i]:
                assignedRangeSize[sortInd[1]] = self.powOfThree[n-i]
            else:
                assignedRangeSize[sortInd[1]] = mpz(gmpy2.ceil(probs[sortInd[1]] * available))
            assignedRangeSize[sortInd[2]] = available - assignedRangeSize[sortInd[1]]
            if seed < L[(i-1)%2] + sum(assignedRangeSize[:1]):
                seq[i-1] = 0
            elif seed < L[(i-1)%2] + sum(assignedRangeSize[:2]):
                seq[i-1] = 1
            elif seed < L[(i-1)%2] + sum(assignedRangeSize[:3]):
                seq[i-1] = 2
            else:
                raise ValueError("Invalid Seed")
            L[i%2] = L[(i-1)%2] + sum(assignedRangeSize[:seq[i-1]])
            U[i%2] = L[(i-1)%2] + sum(assignedRangeSize[:(seq[i-1] + 1)]) - 1
        return U[n%2] - L[n%2] + 1

'''
The DTE model that assumes each sequence has the same probability
'''
class NaiveGenoGuard(BaseGenoGuard):
    def __init__(self, SNPRefList):
        super(NaiveGenoGuard, self).__init__(SNPRefList, mpz(3)**mpz(len(SNPRefList)) - 1)
    '''
    The base class encode method also works for this model, but this re-written implementation is more efficient.
    '''
    def encode(self, seq):
        assert len(seq) == len(self.SNPRefList)
        return gmpy2.mpz(''.join(map(str, seq)), 3)
    def decode(self, seed):
        subSeq = map(int, gmpy2.digits(seed, 3))
        return [0]*(len(self.SNPRefList) - len(subSeq)) + subSeq

'''
The DTE model that uses public linkage disequilibrium and allele frequency dat for building a first-order Markov chain
'''
class PubLDGenoGuard(BaseGenoGuard):
    def __init__(self, SNPRefList, bitLen, AF, LD):
        super(PubLDGenoGuard, self).__init__(SNPRefList, mpz(2)**mpz(bitLen) - 1)
        self.AF = AF
        self.LD = LD
    def condProb(self, prefixSeq):
        prefixLen = len(prefixSeq)
        assert prefixLen >= 1
        SNVLDTuple = self.LD[prefixLen - 1][0:2]
        probs = np.ones((3)) / 3
        # When there is no LD for the current SNV
        if len(SNVLDTuple) == 0:
            probs[0] = self.AF[prefixLen - 1][0] ** 2
            probs[1] = 2 * self.AF[prefixLen - 1][0] * self.AF[prefixLen - 1][1]
            probs[2] = 1 - probs[0] - probs[1]
        else:
            jointMat = self.pairwiseJoint(prefixLen - 1, SNVLDTuple[0], SNVLDTuple[1])
            probs = jointMat[prefixSeq[SNVLDTuple[0]], :]
            probs = probs / sum(probs)
        probs[np.nonzero(probs < self.PREC)] = self.PREC  # To avoid having an extremely low conditional probability, which does not make much sense for real sequences
        probs = probs / sum(probs)    
        return probs
    '''
    Calculate the joint probability of two SNPs, namely, Pr[SNV1, SNV2]
    Actually, SNV1 is the *index* of the current SNV, SNV2 is the *index* of a previous SNV.
    D is the D_value between these two SNVs.
    '''
    def pairwiseJoint(self, SNV1, SNV2, D):
        p_AB = self.AF[SNV2][0] * self.AF[SNV1][0] + D
        p_Ab = self.AF[SNV2][0] * self.AF[SNV1][1] - D
        p_aB = self.AF[SNV2][1] * self.AF[SNV1][0] - D
        p_ab = self.AF[SNV2][1] * self.AF[SNV1][1] + D
        jointMat = np.zeros((3,3))
        jointMat[0, 0] = p_AB * p_AB
        jointMat[0, 1] = 2 * p_AB * p_Ab
        jointMat[0, 2] = p_Ab * p_Ab
        jointMat[1, 0] = 2 * p_AB * p_aB
        jointMat[1, 1] = 2 * p_AB * p_ab + 2 * p_Ab * p_aB
        jointMat[1, 2] = 2 * p_Ab * p_ab
        jointMat[2, 0] = p_aB * p_aB
        jointMat[2, 1] = 2 * p_aB * p_ab
        jointMat[2, 2] = p_ab * p_ab
        jointMat[np.nonzero(jointMat < 0)] = 1e-20
        jointMat = jointMat / sum(sum(jointMat))
        return jointMat

'''
The DTE model that builds k-th order Markov chain based on real genotype datasets
'''
class DirectCondProbGenoGuard(BaseGenoGuard):
    def __init__(self, SNPRefList, bitLen, directCondProbs, order):
        super(DirectCondProbGenoGuard, self).__init__(SNPRefList, mpz(2)**mpz(bitLen) - 1)
        self.directCondProbs = directCondProbs
        self.order = order
    def condProb(self, prefixSeq):
        prefixLen = len(prefixSeq)
        assert prefixLen >= 1
        idx = 0
        for i in range(min(self.order, prefixLen-1), 0, -1):
            idx += prefixSeq[prefixLen - i - 1]
            idx = idx * 3
        probs = np.array(self.directCondProbs[len(prefixSeq) - 1][idx:(idx+3)])
        probs[np.nonzero(probs < self.PREC)] = self.PREC  # To avoid having an extremely low conditional probability, which does not make much sense for real sequences
        probs = probs / sum(probs)
        return probs
        
'''
The DTE model that builds a hidden Markov model based on recombination.
It is actually the implementation of the paper:
J. Marchini, B. Howie, S. Myers, G. McVean, and P. Donnelly, "A new multipoint method for genome-wide association studies by imputation of genotypes", Nature genetics, vol. 39, pp. 906-913, 2007.
'''
class RecombGenoGuard(BaseGenoGuard):
    Ne = 11418   # Effective population size
    def __init__(self, SNPRefList, bitLen, haplotype, geneticDist):
        super(RecombGenoGuard, self).__init__(SNPRefList, mpz(2)**mpz(bitLen) - 1)
        self.bitLen = bitLen
        self.haplotype = np.array(haplotype)
        self.N = np.shape(haplotype)[1]
        self.geneticDist = geneticDist
        self.theta = self.computeTheta(self.N)
        self.mutateMatrix = self.computeMutateMatrix(self.theta, self.N)
        self.alpha = np.zeros((3, self.N, self.N))
    
    '''
    The mutation rate
    '''
    def computeTheta(self, N):
        theta = 0
        for i in range(1, N):
            theta += 1.0 / i
        theta = 1.0 / theta
        theta *= 0.1  # Fine-tuning the mutation rate to fit in real dataset
        return theta
    
    def computeMutateMatrix(self, theta, N):
        r = theta * 0.5 / (theta + N)
        return np.array([[(1-r)*(1-r), 2*r*(1-r), r*r],
                [r*(1-r), r*r + (1-r)*(1-r), r*(1-r)],
                [r*r, 2*r*(1-r), (1-r)*(1-r)]])
        
    def condProb(self, prefixSeq):
        prefixLen = len(prefixSeq)
        prefixProb = 1  # Note that this prefixProb is the probability of the prefix, not including the current snp (unlike prefixLen)
        if prefixLen == 1:
            m1 = np.repeat(np.array([self.haplotype[0, :]]), self.N, 0)
            m2 = np.repeat(np.array([self.haplotype[0, :]]).transpose(), self.N, 1)
            rowIdx = np.array(m1+m2, dtype=np.intp)
            for columnIdx in range(3):
                self.alpha[columnIdx, :, :] = self.mutateMatrix[rowIdx, columnIdx] * 1.0 / (self.N*self.N)
        else:
            # suppose x is a state tuple (x_1, x_2)
            # alpha_{j+1} (x) = r_{j+1} (x) * ( p*p*alpha_j(x) + p*q*sum_{y_1 = x_1}(alpha_j(y)) + p*q*sum_{y_2 = x_2}(alpha_j(y) + q*q*sum_{y}(alpha_j(y)))
#             p = math.exp(-4*self.Ne*self.geneticDist[prefixLen-1] * 1.0 / self.N)
            if self.geneticDist[prefixLen-1] == 0.0:
                self.geneticDist[prefixLen-1] = 1e-8
            p = math.exp(-self.biasCorrection(4*self.Ne*self.geneticDist[prefixLen-1], int(self.SNPRefList[prefixLen-1]) - int(self.SNPRefList[prefixLen - 2])) * 1.0 / self.N)
#             p = math.exp(-4*self.Ne*self.geneticDist[prefixLen-1] * 1.0 / self.N)
            q = (1-p) / self.N
            alpha_j = self.alpha[prefixSeq[prefixLen - 2], :, :]
            # term 1
            term_1 = p*p*alpha_j
            # term 2 and term 3
            sum_vec1 = np.array([np.sum(alpha_j, axis=0)])
            sum_vec2 = np.array([np.sum(alpha_j, axis=1)])
            sum_mat1 = np.repeat(sum_vec1, self.N, 0)
            sum_mat2 = np.repeat(sum_vec2.transpose(), self.N, 1)
            term_2_3 = p*q*(sum_mat1 + sum_mat2)
            # term 4
            prefixProb = np.sum(sum_vec1)
            term_4 = q*q*prefixProb
            
            # alpha_{j+1}
            term = term_1 + term_2_3 + term_4
            m1 = np.repeat(np.array([self.haplotype[prefixLen-1, :]]), self.N, 0)
            m2 = np.repeat(np.array([self.haplotype[prefixLen-1, :]]).transpose(), self.N, 1)
            rowIdx = np.array(m1+m2, dtype=np.intp)
            for columnIdx in range(3):
                self.alpha[columnIdx, :, :] = self.mutateMatrix[rowIdx, columnIdx] * term
            # Renormalize alpha, to avoid float number precision problem when the sequence gets longer and longer, and alpha gets smaller and smaller
            # Since we only need the conditional probabilities, but not the probability of the prefix sequence, it doesn't matter if we scale alpha by a constant.
            # Here we scale it by dividing prefixProb so that alpha sums to 1. But other scale methods also work, only if we maintain the relative ratio between different alpha elements
            self.alpha = self.alpha / prefixProb
        probs = np.sum(np.sum(self.alpha, axis=2), axis=1)
        # Try weighting this with a uniform distribution to avoid overfit
        weight = 1.0
        return weight*probs + (1-weight)*(np.array([0.25, 0.5, 0.25]))
    '''
    Bias correction.
    '''
    def biasCorrection(self, rate, phyDist):
        rho = rate / phyDist
        return rate * 0.0015
            
'''
Load allele frequencies and linkage disequilibrium

AFFileName is the allele frequency file. It has one SNV in each row, with the format: position    major_allele_base    major_allele_frequency    minor_allele_base    minor_allele_frequency
For example, "14795579    G    0.812    T    0.188"

LDFileName is the linkage disequilibrium file. It has one pair of linkage disequilibrium in each row, with the format: current_SNV_position    index_of_the_previous_SNV    D_value
For example, "14805814    0    -0.0326711264643", means SNV 14805814 has LD with the SNV of index 0 (which is SNV 14795579 in the datafile I gave you), with D_value = -0.0326711264643
If a row only has the current_SNV_position column, it means this SNV doesn't have LD with any previous SNV

I define these two formats just for convenience of my code processing. In the future, if we need to consider multiple LDs, we can just add (2*num_of_LDs) columns after the current_SNV_position column.
For example, it could be extended to "15403514    21    -0.195466617566    20    0.132", meaning SNV 15403514 has two LDs, one with SNV of index 21, the other with SNV of index 20.

Please also note that AFFileName and LDFileName should have the same number of rows, namely, the same number of SNVs.
'''
def loadPriorKnowledge(AFFileName, LDFileName):
    AFFile = open(AFFileName)
    AF = []
    for line in AFFile.readlines():
        attrArray = line.split()
        AF.append([float(attrArray[2]), float(attrArray[4])])
    AFFile.close()
    LDFile = open(LDFileName)
    LD = []
    for line in LDFile.readlines():
        attrArray = line.split()
        LD.append(map(lambda x: int(x[1]) if x[0] % 2 == 0 else float(x[1]), enumerate(attrArray[1:])))
    LDFile.close()
    return AF, LD

def loadDirectCondProbs(probFileName):
    probFile = open(probFileName)
    directCondProbs = [map(float, line.split()) for line in probFile.readlines()]
    probFile.close()
    return directCondProbs

def loadRecombData(haplotypeFileName, recombFileName):
    haplotypeFile = open(haplotypeFileName)
    haplotype = map(lambda u: map(int, u.split()[1:]), haplotypeFile.readlines())
    haplotypeFile.close()
    
    recombFile = open(recombFileName)
    geneticDist = map(lambda u: float(u.split()[1]), recombFile.readlines())
    recombFile.close()
    
    return haplotype, geneticDist

if __name__ == "__main__":
    '''
    !IMPORTANT
    Remember to replace the file names with those in your own laptop.
    '''
#     AF, LD = loadPriorKnowledge("../hapmap/small_allele_freqs_chr22_CEU.txt",
#                                 "../hapmap/small_ld_chr22_CEU.txt")
#     directCondProbs = loadDirectCondProbs("../hapmap/small_condProb2_chr22_CEU_ref.txt")  # 2nd-order Markov chain
    haplotype, geneticDist = loadRecombData('../hapmap/chr22/small_CEU.chr22.hap',
                                            '../hapmap/chr22/small_genetic_map_chr22_combined_b36.txt')
    SNPRefList = []
    dataset = []
    datFile = open("../hapmap/chr22/small_genotypes_chr22_CEU.txt")
    for line in datFile.readlines():
        attrArray = line.split()
        SNPRefList.append(attrArray[0])
        dataset.append(attrArray[1:])
    dataset = np.array(dataset, dtype=int)
    bitLen = int(math.ceil(len(SNPRefList)*4))    # 4 is the storage overhead parameter as introduced in the paper.
    bitLen = (bitLen / 128 + 1) * 128    # Just for AES encryption convenience
     
    test_gg = RecombGenoGuard(SNPRefList, bitLen, map(lambda u: u[:234], haplotype), geneticDist)

    seq = list(dataset[:, 1])
    ciphertext = test_gg.hgEncrypt(seq, 'any_password')
    decrypted_seq = test_gg.hgDecrypt(ciphertext, 'any_password')
    print seq == decrypted_seq
