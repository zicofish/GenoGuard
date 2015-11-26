#!/usr/bin/python
'''
Created on Dec 11, 2014

@author: zhihuang
'''
from genRandomSeq import *


'''
Predict the hidden SNVs with the specified model
'''
def predictHiddenSNVsRecombModel(Input_SNVs, model):
    predictedSNVs = np.array([([0]*len(Input_SNVs[0])) for i in range(len(Input_SNVs))])
    for j in range(0, len(Input_SNVs[0])):
        for i in range(0, len(Input_SNVs)):
            probs = model.condProb(predictedSNVs[:(i+1), j])
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

'''
Load the specified model and use it to predict hidden SNVs
NOTE: In this function, I have only implemented loading the RecombModel. 
TODO: For the other two models, PubLDModel and DirectCondProbModel, the extension is very simple. Please refer to 
the function 'genRandSeq' in the file 'genRandomSeq.py' to understand how to load the other two models.
'''
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
    # Here we only predict for one person. But inputSNVs is two-dimension, so you can always predict for multiple persons together.
    inputSNVs = map(lambda u: u.split()[(patientIdx+1):(patientIdx+2)], inputSNVsFile.readlines()) 
    
    print "Predicting hidden SNVs for patient ", patientIdx, "..."
    predictedSNVs = predictHiddenSNVsRecombModel(inputSNVs, model)
    
    predictedSNVsFile = open(kwargs['predictedSNVsFileName'], 'w')
    predictedSNVsFile.write("\n".join(map(lambda u: " ".join(map(str, u)), predictedSNVs.transpose())))
    predictedSNVsFile.close()
    print "Complete. Result written to ", kwargs['predictedSNVsFileName']
    

if __name__ == '__main__':
    ####  For running on the cluster   ########
#     import sys
#     arg = sys.argv
#     arg_dict = dict(map(lambda u: u.split('='), arg[2:]))
#     current_module = sys.modules[__name__]
#     getattr(current_module, arg[1])(arg_dict)
    ####  cluster end   #########
    
    loadAndPredict({'genotypeFileName':'../hapmap/chr22/small_genotypes_chr22_CEU.txt',
                    'haplotypeFileName':'../hapmap/chr22/small_CEU.chr22.hap',
                    'geneticDistFileName':'../hapmap/chr22/small_genetic_map_chr22_combined_b36.txt',
                    'hiddenSNVsFileName':'../hapmap/chr22/hiddenSNVs/hiddenSNVs_s0_e165_chr22_0.3.txt',
                    'patientIdx':0,
                    'predictedSNVsFileName':'../hapmap/chr22/hiddenSNVs/predictedSNVs_patient0_0.3.txt'})
