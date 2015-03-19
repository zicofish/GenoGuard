#!/usr/bin/python
'''
Created on Dec 11, 2014

@author: zhihuang

The difference between this file and inference.py is: we don't import genRandomSeq.py here
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
        self.a_ = -3.817e-01 + 6.350e-03 * len(SNPRefList) - 3.833e-05 * len(SNPRefList) * len(SNPRefList);
        self.b_ = -1.133e-01 - 2.600e-03 * len(SNPRefList) + 1.333e-05 * len(SNPRefList) * len(SNPRefList);
        print "a_: ", self.a_
        print "b_: ", self.b_
        
    def computeTheta(self, N):
        theta = 0
        for i in range(1, N):
            theta += 1.0 / i
        theta = 1.0 / theta
        theta *= 0.1
        return theta
    
    def computeMutateMatrix(self, theta, N):
        r = theta * 0.5 / (theta + N)
        print r
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


def loadRecombData(haplotypeFileName, recombFileName):
    haplotypeFile = open(haplotypeFileName)
    haplotype = map(lambda u: map(int, u.split()[1:]), haplotypeFile.readlines())
    haplotypeFile.close()
    
    recombFile = open(recombFileName)
    geneticDist = map(lambda u: float(u.split()[1]), recombFile.readlines())
    recombFile.close()
    
    return haplotype, geneticDist


def predictHiddenSNVsRecombModel(Input_SNVs, recombModel):
    predictedSNVs = np.array([([0]*len(Input_SNVs[0])) for i in range(len(Input_SNVs))])
    for j in range(0, len(Input_SNVs[0])):
        for i in range(0, len(Input_SNVs)):
            probs = recombModel.condProb(predictedSNVs[:(i+1), j])
            if Input_SNVs[i][j] != 'H':
                predictedSNVs[i, j] = int(Input_SNVs[i][j])
            else:
                r = random.random()
                if r < probs[0]:
                    predictedSNVs[i, j] = 0
                elif r < probs[0]+probs[1]:
                    predictedSNVs[i, j] = 1
                elif r < probs[0]+probs[1]+probs[2]:
                    predictedSNVs[i, j] = 2
                else:
                    raise ValueError("Invalid probs: %s, %s, %s" % (probs[0], probs[1], probs[2]))
                
    return predictedSNVs

def loadAndPredict(kwargs):
    SNPRefList = []
    dataset = []
    datFile = open(kwargs['genotypeFileName'])
    for line in datFile.readlines():
        attrArray = line.split()
        SNPRefList.append(attrArray[0])
        dataset.append(attrArray[1:])
    dataset = np.array(dataset, dtype=int)
    
    haplotype, geneticDist = loadRecombData(kwargs['haplotypeFileName'], kwargs['geneticDistFileName'])
    model = RecombModel(SNPRefList, map(lambda u: u[:146], haplotype), geneticDist)
    
    inputSNVsFile = open(kwargs['hiddenSNVsFileName'])
    patientIdx = int(kwargs['patientIdx'])
    inputSNVs = map(lambda u: u.split()[patientIdx:(patientIdx+1)], inputSNVsFile.readlines())
    
    predictedSNVs = predictHiddenSNVsRecombModel(inputSNVs, model)
    
    predictedSNVsFile = open(kwargs['predictedSNVsFileName'], 'w')
    predictedSNVsFile.write("\n".join(map(lambda u: " ".join(map(str, u)), predictedSNVs.transpose())))
    predictedSNVsFile.close()
    

if __name__ == '__main__':
    ####  cluster only   ########
    import sys
    arg = sys.argv
    arg_dict = dict(map(lambda u: u.split('='), arg[2:]))
    current_module = sys.modules[__name__]
    getattr(current_module, arg[1])(arg_dict)
    ####  cluster end   #########
#     loadAndPredict({'genotypeFileName':'./hapmap/chr22/big_genotypes_chr22_CEU.txt',
#                     'haplotypeFileName':'./hapmap/chr22/big_CEU.chr22.hap',
#                     'geneticDistFileName':'./hapmap/chr22/big_genetic_map_chr22_combined_b36.txt',
#                     'hiddenSNVsFileName':'./hapmap/chr22/HiddenSNVs10%1.txt',
#                     'patientIdx':0,
#                     'predictedSNVsFileName':'./hapmap/chr22/predictedSNVs10%1.txt'})
#     