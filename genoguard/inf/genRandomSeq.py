'''
Created on Jun 17, 2014

@author: zhihuang
'''
import math
import random
import numpy as np

class BaseModel(object):
    PREC = 10**-4
    def __init__(self, SNPRefList):
        self.SNPRefList = SNPRefList
        
    def gen(self):
        n = len(self.SNPRefList)
        seq = [-1] * n
        for i in range(1, n + 1):
            probs = self.condProb(seq[:i])
            r = random.random()
            if r < probs[0]:
                seq[i-1] = 0
            elif r < probs[0] + probs[1]:
                seq[i-1] = 1
            elif r < 1:
                seq[i-1] = 2
            else:
                raise ValueError("Invalid random float number: %s" % (r))
        return seq  
    
    def condProb(self, prefixSeq):
        return [0, 0, 0]
    

class PubLDModel(BaseModel):
    def __init__(self, SNPRefList, AF, LD):
        super(PubLDModel, self).__init__(SNPRefList)
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
        probs[np.nonzero(probs < self.PREC)] = self.PREC
        probs = probs / sum(probs)    
        return probs
    '''
    Calculate the joint probability of two SNPs, namely, Pr[SNV1, SNV2]
    In our case, varSNP is the index of the current SNV, condSNP is the index of a previous SNV.
    D is the D_value between these two SNVs.
    '''
    def pairwiseJoint(self, SNV1, SNV2, D):
    #     D = LD[condSNP][varSNP]
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
#         jointMat[np.nonzero(jointMat < self.PREC)] = self.PREC
        jointMat = jointMat / sum(sum(jointMat))
        return jointMat

class DirectCondProbModel(BaseModel):
    def __init__(self, SNPRefList, directCondProbs, order):
        super(DirectCondProbModel, self).__init__(SNPRefList)
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
        probs[np.nonzero(probs < self.PREC)] = self.PREC
        probs = probs / sum(probs)
        return probs
        

class RecombModel(BaseModel):
    Ne = 11418   # Effective population size
    def __init__(self, SNPRefList, haplotype, geneticDist):
        super(RecombModel, self).__init__(SNPRefList)
        self.haplotype = np.array(haplotype)
        self.N = np.shape(haplotype)[1]
        self.geneticDist = geneticDist
        self.theta = self.computeTheta(self.N)
        self.mutateMatrix = self.computeMutateMatrix(self.theta, self.N)
        self.alpha = np.zeros((3, self.N, self.N))
        # Bias correction coefficients
        # From condlike.hpp in "http://stephenslab.uchicago.edu/software.html#hotspotter"
        self.a_ = -3.817e-01 + 6.350e-03 * len(SNPRefList) - 3.833e-05 * len(SNPRefList) * len(SNPRefList);
        self.b_ = -1.133e-01 - 2.600e-03 * len(SNPRefList) + 1.333e-05 * len(SNPRefList) * len(SNPRefList);
        
    def computeTheta(self, N):
        theta = 0
        for i in range(1, N):
            theta += 1.0 / i
        theta = 1.0 / theta
        theta *= 0.1
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
    def biasCorrection(self, rate, phyDist):
        rho = rate / phyDist
        return rate * 0.0015
#         return rate*pow (10.0, self.a_ + self.b_ * math.log10 (rho))
            
'''
Load allele frequencies (AF) and linkage disequilibrium (LD) for the PubLDModel

AFFileName is the allele frequency file. It has one SNV in each row, with the format: position    major_allele_base    major_allele_frequency    minor_allele_base    minor_allele_frequency
For example, "14795579    G    0.812    T    0.188"

LDFileName is the linkage disequilibrium file. It has one pair of linkage disequilibrium in each row, with the format: current_SNV_position    index_of_the_previous_SNV    D_value
For example, "14805814    0    -0.0326711264643", means SNV 14805814 has LD with the SNV of index 0 (which is SNV 14795579 in the datafile I gave you), with D_value = -0.0326711264643
If a row only has the current_SNV_position column, it means this SNV doesn't have LD with any previous SNV

I define these two formats just for convenience of my code processing. In the future, if we need to consider multiple LDs, we can just add (2*num_of_LDs) columns after the current_SNV_position column.
For example, it could be extended to "15403514    21    -0.195466617566    20    0.132", meaning SNV 15403514 has two LDs, one with SNV of index 21, the other with SNV of index 20.

Please also note that AFFileName and LDFileName should have the same number of rows, namely, the same number of SNVs.
'''
def loadAFAndLD(AFFileName, LDFileName):
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
Load k-th order markov chain for the DirectCondProbModel
'''
def loadDirectCondProbs(probFileName):
    probFile = open(probFileName)
    directCondProbs = [map(float, line.split()) for line in probFile.readlines()]
    probFile.close()
    return directCondProbs

'''
Load haplotype and genetic distance for the RecombModel
'''
def loadRecombData(haplotypeFileName, recombFileName):
    haplotypeFile = open(haplotypeFileName)
    haplotype = map(lambda u: map(int, u.split()[1:]), haplotypeFile.readlines())
    haplotypeFile.close()
    
    recombFile = open(recombFileName)
    geneticDist = map(lambda u: float(u.split()[1]), recombFile.readlines())
    recombFile.close()
    
    return haplotype, geneticDist

'''
Generate random samples under the specified model
'''
def genRandSeq(kwargs):
    SNPRefList = []
    dataset = []
    datFile = open(kwargs['genotypeFileName'])
    for line in datFile.readlines():
        attrArray = line.split()
        SNPRefList.append(attrArray[0])
        dataset.append(attrArray[1:])
    dataset = np.array(dataset, dtype=int)
    
    if kwargs['modelName'] == 'PubLDModel':
        AF, LD = loadAFAndLD(kwargs['AFFileName'], kwargs['LDFileName'])
        model = PubLDModel(SNPRefList, AF, LD)
    elif kwargs['modelName'] == 'DirectCondProbModel':
        directCondProbs = loadDirectCondProbs(kwargs['probFileName'])
        model = DirectCondProbModel(SNPRefList, directCondProbs, int(kwargs['order']))
    elif kwargs['modelName'] == 'RecombModel':
        haplotype, geneticDist = loadRecombData(kwargs['haplotypeFileName'],
                                                kwargs['geneticDistFileName'])
        model = RecombModel(SNPRefList, map(lambda u: u[:200], haplotype), geneticDist)
    outFile = open(kwargs['outFileName'], 'a')
    for i in range(int(kwargs['sampleSize'])):
        print "Generating sample sequence ", i, "..."
        seq = model.gen()
        outFile.write(' '.join(map(str, seq)) + '\n')
    outFile.close()

if __name__ == "__main__":
    '''
    !IMPORTANT
    Remember to replace the file names with those in your own laptop.
    '''
#     genRandSeq({'genotypeFileName':"../hapmap/chr22/small_genotypes_chr22_CEU.txt",
#                 'modelName':'DirectCondProbModel',
#                 'probFileName':'../hapmap/chr22/small_condProb0_chr22_CEU_ref.txt',
#                 'order':0,
#                 'outFileName':'../hapmap/chr22/small_directcond_random_chr22_CEU.txt',
#                 'sampleSize':100})
    genRandSeq({'genotypeFileName':"../hapmap/chr22/small_genotypes_chr22_CEU.txt",
                'modelName':'RecombModel',
                'haplotypeFileName':'../hapmap/chr22/small_CEU.chr22.hap',
                'geneticDistFileName':'../hapmap/chr22/small_genetic_map_chr22_combined_b36.txt',
                'outFileName':'../hapmap/chr22/small_recomb_random_chr22_CEU.txt',
                'sampleSize':10})