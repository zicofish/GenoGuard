'''
Created on Nov 18, 2014

@author: zhihuang
'''

import numpy as np

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
    jointMat[np.nonzero(jointMat < 10**-8)] = 10**-8
    jointMat = jointMat / sum(sum(jointMat))
    return jointMat

def ldModel(AF, LD):
    cond_ld = []
    for i in range(len(AF)):
        if LD[i] == []:
            p_0 = AF[i][0] * AF[i][0]
            p_1 = AF[i][0] * AF[i][1]
            p_2 = AF[i][1] * AF[i][1]
            cond_ld.append([p_0, p_1, p_2, p_0, p_1, p_2, p_0, p_1, p_2])
        else:
            joinMat = pairwiseJoint(i, LD[i][0], AF, LD[i][1])
            prefixProb = np.sum(joinMat, axis=1)
            ext_prefixProb = np.repeat(np.array([prefixProb]).transpose(), 3, axis=1)
            cond_probs = joinMat / ext_prefixProb
            cond_ld.append(list(cond_probs.flatten()))
    return cond_ld

if __name__ == "__main__":
    '''
    !IMPORTANT
    Remember to replace the file names with those in your own laptop.
    '''
    AF, LD = loadPriorKnowledge("./hapmap/chr8/big_allele_freqs_chr8_CEU.txt",
                                "./hapmap/chr8/big_ld_chr8_CEU.txt")
    model = ldModel(AF, LD)
    tmpfile = open("./hapmap/chr8/LDModel_chr8_CEU.txt", 'w')
    tmpfile.write('\n'.join(map(lambda u: '\t'.join(map(str, u)), model)))
    tmpfile.close()