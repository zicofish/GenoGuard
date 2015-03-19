#!/usr/bin/python
'''
Created on Dec 11, 2014

@author: zhihuang
'''
from genRandomSeq import *
import numpy as np

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