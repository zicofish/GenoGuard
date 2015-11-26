'''
Created on Jun 17, 2014

@author: zhihuang

@deprecated: This file is an initial product of the solution. It has been replaced by GenoGuard.py.
'''
import math
import random
import numpy as np
from gmpy2 import mpz
import struct
from Crypto.Hash import HMAC
from Crypto.Cipher import AES
from Crypto import Random
import time

PREC = 10**9
powOfThree = []


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

'''
The coding function.
Please do not modify it. Contact me first if you need to change something.
'''
def encode(seq, SNPRefList, AF, LD, bitLen = 1000):
    assert len(seq) == len(SNPRefList)
    assert 2**bitLen >= 3**len(seq)
    
    L = [-1, -1]; U = [-1, -1]
    L[0] = mpz(0); U[0] = mpz(2**bitLen - 1)
    n = len(seq)
    for i in range(1, n + 1):
        probs = condProb(seq[:i], SNPRefList, AF, LD)
        sortInd = np.argsort(probs)
        assignedRangeSize = [0] * 3
        available = U[(i-1)%2] - L[(i-1)%2] + 1
        if int(probs[sortInd[0]] * PREC) \
        * available / PREC < powOfThree[n-i]:
            assignedRangeSize[sortInd[0]] = powOfThree[n-i]
        else:
            assignedRangeSize[sortInd[0]] = (int(probs[sortInd[0]] * PREC) \
            * available + PREC -1 ) / PREC
        available -= assignedRangeSize[sortInd[0]]
        probs[sortInd[0]] = 0
        probs = probs / sum(probs)  #re-normalize the remaining two probabilities
        if int(probs[sortInd[1]] * PREC) \
        * available / PREC < powOfThree[n-i]:
            assignedRangeSize[sortInd[1]] = powOfThree[n-i]
        else:
            assignedRangeSize[sortInd[1]] = (int(probs[sortInd[1]] * PREC) \
            * available + PREC -1 ) / PREC
        assignedRangeSize[sortInd[2]] = available - assignedRangeSize[sortInd[1]]
        L[i%2] = L[(i-1)%2] + sum(assignedRangeSize[:seq[i-1]])
        U[i%2] = L[(i-1)%2] + sum(assignedRangeSize[:(seq[i-1]+1)]) - 1
    seed = random.randint(L[n%2], U[n%2])
    return seed

'''
The decoding function. Almost symmetric to coding.
'''
def decode(seed, SNPRefList, AF, LD, bitLen = 1000):
    L = [-1, -1]; U = [-1, -1]
    L[0] = mpz(0); U[0] = mpz(2**bitLen - 1)
    n = len(SNPRefList)
    seq = [-1] * n
    for i in range(1, n + 1):
        probs = condProb(seq[:i], SNPRefList, AF, LD)
        sortInd = np.argsort(probs)
        assignedRangeSize = [0] * 3
        available = U[(i-1)%2] - L[(i-1)%2] + 1
        if int(probs[sortInd[0]] * PREC) \
        * available / PREC < powOfThree[n-i]:
            assignedRangeSize[sortInd[0]] = powOfThree[n-i]
        else:
            assignedRangeSize[sortInd[0]] = (int(probs[sortInd[0]] * PREC) \
            * available + PREC -1 ) / PREC
        available -= assignedRangeSize[sortInd[0]]
        probs[sortInd[0]] = 0
        probs = probs / sum(probs)  #re-normalize the remaining two probabilities
        if int(probs[sortInd[1]] * PREC) \
        * available / PREC < powOfThree[n-i]:
            assignedRangeSize[sortInd[1]] = powOfThree[n-i]
        else:
            assignedRangeSize[sortInd[1]] = (int(probs[sortInd[1]] * PREC) \
            * available + PREC -1 ) / PREC
        assignedRangeSize[sortInd[2]] = available - assignedRangeSize[sortInd[1]]
        if seed < L[(i-1)%2] + sum(assignedRangeSize[:1]):
            seq[i-1] = 0
        elif seed < L[(i-1)%2] + sum(assignedRangeSize[:2]):
            seq[i-1] = 1
        elif seed < L[(i-1)%2] + sum(assignedRangeSize[:3]):
            seq[i-1] = 2
        else:
            raise ValueError("Invalid Seed: " + seed)
        L[i%2] = L[(i-1)%2] + sum(assignedRangeSize[:seq[i-1]])
        U[i%2] = L[(i-1)%2] + sum(assignedRangeSize[:(seq[i-1] + 1)]) - 1
    return seq

'''
Based on kth-order markov chain, calculate the conditional probability of current SNV, given k previous SNV.
Here, current SNV is prefixSeq[ prefixLen - 1 ], namely, the last element of prefixSeq. And prefixSeq contains all SNVs before the current SNV, from the first one to the current one.
SNPRefList is the list of all SNV positions, from the first one to the last one.
AF and LD is what we have loaded in loadPriorKnowledge.
'''
def condProb(prefixSeq, SNPRefList, AF, LD):
    condSNPs = heuristicCondSNPs(prefixSeq, SNPRefList, LD)
    probs = np.ones((3)) / 3
    prefixLen = len(prefixSeq)
    condSNPs=[]
    # When there is no LD for the current SNV
    if len(condSNPs) == 0:
        probs[0] = math.floor(AF[prefixLen - 1][0] ** 2 \
                              * PREC) / PREC
        probs[1] = math.floor(2 * AF[prefixLen - 1][0] \
                              * AF[prefixLen - 1][1] \
                              * PREC) / PREC
        probs[2] = 1 - probs[0] - probs[1]
        
    # When there is one LD for the current SNV
    elif len(condSNPs) == 2:
        jointMat = pairwiseJoint(prefixLen - 1, condSNPs[0], AF, condSNPs[1])
        probs = jointMat[prefixSeq[condSNPs[0]], :]
        probs = probs / sum(probs)
        probs[0] = math.floor(probs[0] * PREC) / PREC
        probs[1] = math.floor(probs[1] * PREC) / PREC
        probs[2] = 1 - probs[0] - probs[1]
        
    # When there is more than one LD...
    
    return probs

'''
Return the previous SNVs which we want to take into account when dealing with the kth-order markov chain.
As what we did in our paper, we only return the previous one. But we can adapt this function to take more into account.
'''
def heuristicCondSNPs(prefixSeq, SNPRefList, LD):
    """
    Only consider the nearest SNP
    """
    assert len(prefixSeq) >= 1
    return LD[len(prefixSeq) - 1][0:2]
    
'''
Calculate the joint probability of two SNPs, namely, Pr[SNV1, SNV2]
In our case, varSNP is the index of the current SNV, condSNP is the index of a previous SNV.
D is the D_value between these two SNVs.
'''
def pairwiseJoint(varSNP, condSNP, AF, D):
#     D = LD[condSNP][varSNP]
    p_AB = AF[condSNP][0] * AF[varSNP][0] + D
    p_Ab = AF[condSNP][0] * AF[varSNP][1] - D
    p_aB = AF[condSNP][1] * AF[varSNP][0] - D
    p_ab = AF[condSNP][1] * AF[varSNP][1] + D
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
#     jointMat[np.nonzero(jointMat < 0)] = 0
    jointMat[np.nonzero(jointMat < 10**-4)] = 10**-4
    jointMat = jointMat / sum(sum(jointMat))
    return jointMat

'''
Pack a big integer into binary dat.
This is only used for AES block cipher. You don't need this.
'''
def int_to_packed(val, bitLen):
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
This is only used for AES block cipher. You don't need this.
'''
def packed_to_int(packed, bitLen):
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
You don't need this.
'''
def pbeEncrypt(msg, password):
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
You don't need this.
'''
def pbeDecrypt(ext_c, password):
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
You don't need this.
'''
def hgEncrypt(seq, SNPRefList, AF, LD, bitLen, password):
    seed = encode(seq, SNPRefList, AF, LD, bitLen)
    msg = int_to_packed(seed, bitLen)
    ext_c = pbeEncrypt(msg, password)
    return ext_c

'''
Honey decryption including password-based decryption and decoding.
You don't need this.
'''
def hgDecrypt(ext_c, SNPRefList, AF, LD, bitLen, password):
    msg = pbeDecrypt(ext_c, password)
    seed = packed_to_int(msg, bitLen)
    seq = decode(seed, SNPRefList, AF, LD, bitLen)
    return seq

'''
You don't need this.
'''    
def genRandomSamples(ancestry, bitLen, SNPRefList):
    test_AF, test_LD = loadPriorKnowledge(ancestry + "_freq.txt", "")
    experiments = []
    for i in range(100):
        randomSeed = random.randint(0, 2**bitLen-1)
        randomSeq = decode(randomSeed, SNPRefList, test_AF, test_LD, bitLen)
        experiments.append(['EXP'+str(i), 'TEST'+ancestry] + map(str, randomSeq))
    expFile = open("exp" + ancestry + ".txt", 'w')
    expFile.write('\n'.join(['\t'.join(randSample) for randSample in experiments]))
    expFile.close()

if __name__ == "__main__":
    '''
    !IMPORTANT
    Remember to replace the file names with those in your own laptop.
    '''
    AF, LD = loadPriorKnowledge("C:/Users/zhihuang/Desktop/dat/hapmap/hg_allele_freqs_chr22_CEU.txt",
                                "C:/Users/zhihuang/Desktop/dat/hapmap/hg_ld_chr22_CEU.txt")
    SNPRefList = []
    dataset = []
    datFile = open("C:/Users/zhihuang/Desktop/dat/hapmap/new_genotypes_chr22_CEU.txt")
    for line in datFile.readlines():
        attrArray = line.split()
        SNPRefList.append(attrArray[0])
        dataset.append(attrArray[1:])
    dataset = np.array(dataset, dtype=int)
    bitLen = int(math.ceil(len(SNPRefList)*3.5))    # 3.5 is the storage overhead parameter as introduced in the paper.
    bitLen = (bitLen / 128 + 1) * 128    # Just for AES encryption convenience
    for i in range(len(SNPRefList)):    # For efficiency, pre-compute large powers of 3, as they will be used heavily in encoding and decoding.
        powOfThree.append(mpz(3) ** mpz(i))
        
    seq = list(dataset[:, 0])   # The first sequence in the dataset
    seed = encode(seq, SNPRefList, AF, LD, bitLen)
    print seed  # You may not see it in the console, because it is a huge integer. If you indeed want to see this big integer, please play with a shorter sequence.
    m_seq = decode(seed, SNPRefList, AF, LD, bitLen)
    print seq == m_seq
    