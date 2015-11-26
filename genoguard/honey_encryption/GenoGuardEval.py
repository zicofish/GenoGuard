#!/usr/bin/python
'''
Created on Jul 16, 2014

@author: zhihuang
'''
from GenoGuard import *

def plotResult(datasetProbFileName, randomProbFileName):
    datasetProbFile = open(datasetProbFileName)
    datasetProbs = map(float, datasetProbFile.readlines())
    datasetProbFile.close()
    
    randomProbFile = open(randomProbFileName)
    randomProbs = map(float, randomProbFile.readlines())
    randomProbFile.close()
    
    randomMin = min(randomProbs)
    randomMax = max(randomProbs)
    print (sum(np.array(datasetProbs) < randomMin) + sum(np.array(datasetProbs) > randomMax))*1.0 / 165
    
#     colors = ['blue'] * len(naiveResult)
#     colors[orgInd] = 'red'
#     sizes = [50] * len(naiveResult)
#     sizes[orgInd] = 150
    
    import matplotlib.pyplot as plt
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    h1 = ax.scatter(range(165), datasetProbs, marker='<', s = 100, c='red')
    h2 = ax.scatter([random.randint(0, 164) for i in range(len(randomProbs))], randomProbs, marker='o', s = 50, c = 'blue')
    plt.show()
    
#     fig = plt.figure()
#     ax1 = fig.add_subplot(1, 2, 1)
#     h1 = ax1.scatter(pwds, naiveResult, marker='o', s = 50, c='blue')
#     h2 = ax1.scatter([correctPwd], [correctValue], marker='<', s = 500, c = 'red')
#     ax1.set_xlabel('Password', fontsize=16)
#     ax1.set_ylabel('Logarithm of Interval Size', fontsize=16)
#     ax1.legend((h1, h2),
#            ('Wrong Sequence', 'Correct Sequence'),
#            bbox_to_anchor=(0.50, 1.02, 1., .102),
#            scatterpoints=1,
#            loc='lower right',
#            ncol=2,
#            fontsize=12)
#     
#     ax2 = fig.add_subplot(1, 2, 2)
#     ax2.scatter(pwds, honeyResult, marker='o', s = 50, c='blue')
#     ax2.scatter([correctPwd], [correctValue], marker='<', s = 500, c = 'red')
#     ax2.set_xlabel('Password', fontsize=16)
#     ax2.set_ylabel('Logarithm of Interval Size', fontsize=16)
#     

def generateModelSample(kwargs):
    SNPRefList = []
    dataset = []
    datFile = open(kwargs['genotypeFileName'])
    for line in datFile.readlines():
        attrArray = line.split()
        SNPRefList.append(attrArray[0])
        dataset.append(attrArray[1:])
    dataset = np.array(dataset, dtype=int)
    bitLen = int(math.ceil(len(SNPRefList)*2.5))    # 3.5 is the storage overhead parameter as introduced in the paper.
    bitLen = (bitLen / 128 + 1) * 128    # Just for AES encryption convenience
    
    if kwargs['modelName'] == 'PubLDGenoGuard':
        test_AF, test_LD = loadPriorKnowledge(kwargs['AFFileName'], kwargs['LDFileName'])
        test_gg = PubLDGenoGuard(SNPRefList, bitLen, test_AF, test_LD)
    elif kwargs['modelName'] == 'DirectCondProbGenoGuard':
        directCondProbs = loadDirectCondProbs(kwargs['probFileName'])
        test_gg = DirectCondProbGenoGuard(SNPRefList, bitLen, directCondProbs, int(kwargs['order']))
    outFile = open(kwargs['outFileName'], 'w')
    randomState = gmpy2.random_state(random.randint(0, 2**30))
    for i in range(int(kwargs['sampleSize'])):
        randomSeed = gmpy2.mpz_urandomb(randomState, bitLen)
        randomIntervalSize = test_gg.seedIntervalSize(randomSeed)
        outFile.write(str(gmpy2.log2(randomIntervalSize)) + '\n')
    outFile.close()
    
def computeDatasetProbs(kwargs):
    SNPRefList = []
    dataset = []
    datFile = open(kwargs['genotypeFileName'])
    for line in datFile.readlines():
        attrArray = line.split()
        SNPRefList.append(attrArray[0])
        dataset.append(attrArray[1:])
    dataset = np.array(dataset, dtype=int)
    bitLen = int(math.ceil(len(SNPRefList)*2.5))    # 3.5 is the storage overhead parameter as introduced in the paper.
    bitLen = (bitLen / 128 + 1) * 128    # Just for AES encryption convenience
    
    if kwargs['modelName'] == 'PubLDGenoGuard':
        test_AF, test_LD = loadPriorKnowledge(kwargs['AFFileName'], kwargs['LDFileName'])
        test_gg = PubLDGenoGuard(SNPRefList, bitLen, test_AF, test_LD)
    elif kwargs['modelName'] == 'DirectCondProbGenoGuard':
        directCondProbs = loadDirectCondProbs(kwargs['probFileName'])
        test_gg = DirectCondProbGenoGuard(SNPRefList, bitLen, directCondProbs, int(kwargs['order']))
        
    outFile = open(kwargs['outFileName'], 'w')
    for i in range(int(kwargs['idxBegin']), int(kwargs['idxEnd'])):
        intervalSize = test_gg.seqIntervalSize(list(dataset[:, i]))
        outFile.write(str(gmpy2.log2(intervalSize)) + '\n')
    outFile.close()
        
if __name__ == "__main__":
#     import sys
#     arg = sys.argv
#     arg_dict = dict(map(lambda u: u.split('='), arg[2:]))
#     current_module = sys.modules[__name__]
#     getattr(current_module, arg[1])(arg_dict)

###############################   small    ###########################################
    plotResult('./tmp/publd_dataset0.txt', './tmp/publd_sample0.txt')
    plotResult('./tmp/directcondprob_plain_dataset1.txt', './tmp/directcondprob_plain_sample1.txt')
    plotResult('./tmp/directcondprob_plain_dataset2.txt', './tmp/directcondprob_plain_sample2.txt')
#     computeDatasetProbs({'modelName':'PubLDGenoGuard',
#                          'genotypeFileName':'./hapmap/small_genotypes_chr22_CEU.txt',
#                          'AFFileName':'./hapmap/small_allele_freqs_chr22_CEU.txt',
#                          'LDFileName':'./hapmap/small_ld_chr22_CEU.txt',
#                          'outFileName':'./tmp/publd_dataset0.txt',
#                          'idxBegin':0,
#                          'idxEnd':165})
#     generateModelSample({'modelName':'PubLDGenoGuard',
#                         'genotypeFileName':'./hapmap/small_genotypes_chr22_CEU.txt',
#                         'AFFileName':'./hapmap/small_allele_freqs_chr22_CEU.txt',
#                         'LDFileName':'./hapmap/small_ld_chr22_CEU.txt',
#                         'outFileName':'./tmp/publd_sample0.txt',
#                         'sampleSize':500})
#     computeDatasetProbs({'modelName':'DirectCondProbGenoGuard',
#                          'genotypeFileName':'./hapmap/small_genotypes_chr22_CEU.txt',
#                          'probFileName':'./hapmap/condProb0_chr22_CEU_ref.txt',
#                          'outFileName':'./tmp/directcondprob_plain_dataset0.txt',
#                          'order':0,
#                          'idxBegin':0,
#                          'idxEnd':165})
#     generateModelSample({'modelName':'DirectCondProbGenoGuard',
#                         'genotypeFileName':'./hapmap/small_genotypes_chr22_CEU.txt',
#                         'probFileName':'./hapmap/condProb0_chr22_CEU_ref.txt',
#                         'outFileName':'./tmp/directcondprob_plain_sample0.txt',
#                         'order':0,
#                         'sampleSize':500})
#     computeDatasetProbs({'modelName':'DirectCondProbGenoGuard',
#                          'genotypeFileName':'./hapmap/small_genotypes_chr22_CEU.txt',
#                          'probFileName':'./hapmap/condProb1_chr22_CEU_ref.txt',
#                          'outFileName':'./tmp/directcondprob_plain_dataset1.txt',
#                          'order':1,
#                          'idxBegin':0,
#                          'idxEnd':165})
#     generateModelSample({'modelName':'DirectCondProbGenoGuard',
#                         'genotypeFileName':'./hapmap/small_genotypes_chr22_CEU.txt',
#                         'probFileName':'./hapmap/condProb1_chr22_CEU_ref.txt',
#                         'outFileName':'./tmp/directcondprob_plain_sample1.txt',
#                         'order':1,
#                         'sampleSize':500})
#     computeDatasetProbs({'modelName':'DirectCondProbGenoGuard',
#                          'genotypeFileName':'./hapmap/small_genotypes_chr22_CEU.txt',
#                          'probFileName':'./hapmap/condProb2_chr22_CEU_ref.txt',
#                          'outFileName':'./tmp/directcondprob_plain_dataset2.txt',
#                          'order':2,
#                          'idxBegin':0,
#                          'idxEnd':165})
#     generateModelSample({'modelName':'DirectCondProbGenoGuard',
#                         'genotypeFileName':'./hapmap/small_genotypes_chr22_CEU.txt',
#                         'probFileName':'./hapmap/condProb2_chr22_CEU_ref.txt',
#                         'outFileName':'./tmp/directcondprob_plain_sample2.txt',
#                         'order':2,
#                         'sampleSize':500})
#     computeDatasetProbs({'modelName':'DirectCondProbGenoGuard',
#                          'genotypeFileName':'./hapmap/small_genotypes_chr22_CEU.txt',
#                          'probFileName':'./hapmap/imm2_new_chr22_CEU_ref.txt',
#                          'outFileName':'./tmp/directcondprob_dataset2.txt',
#                          'order':2,
#                          'idxBegin':0,
#                          'idxEnd':10})
#     generateModelSample({'modelName':'DirectCondProbGenoGuard',
#                         'genotypeFileName':'./hapmap/small_genotypes_chr22_CEU.txt',
#                         'probFileName':'./hapmap/imm2_new_chr22_CEU_ref.txt',
#                         'outFileName':'./tmp/directcondprob_sample2.txt',
#                         'order':2,
#                         'sampleSize':10})
#     computeDatasetProbs({'modelName':'DirectCondProbGenoGuard',
#                          'genotypeFileName':'./hapmap/small_genotypes_chr22_CEU.txt',
#                          'probFileName':'./hapmap/imm1_new_chr22_CEU_ref.txt',
#                          'outFileName':'./tmp/directcondprob_dataset1.txt',
#                          'order':1,
#                          'idxBegin':0,
#                          'idxEnd':10})
#     generateModelSample({'modelName':'DirectCondProbGenoGuard',
#                         'genotypeFileName':'./hapmap/small_genotypes_chr22_CEU.txt',
#                         'probFileName':'./hapmap/imm1_new_chr22_CEU_ref.txt',
#                         'outFileName':'./tmp/directcondprob_sample1.txt',
#                         'order':1,
#                         'sampleSize':10})
#     computeDatasetProbs({'modelName':'DirectCondProbGenoGuard',
#                          'genotypeFileName':'./hapmap/small_genotypes_chr22_CEU.txt',
#                          'probFileName':'./hapmap/imm0_new_chr22_CEU_ref.txt',
#                          'outFileName':'./tmp/directcondprob_dataset0.txt',
#                          'order':0,
#                          'idxBegin':0,
#                          'idxEnd':10})
#     generateModelSample({'modelName':'DirectCondProbGenoGuard',
#                         'genotypeFileName':'./hapmap/small_genotypes_chr22_CEU.txt',
#                         'probFileName':'./hapmap/imm0_new_chr22_CEU_ref.txt',
#                         'outFileName':'./tmp/directcondprob_sample0.txt',
#                         'order':0,
#                         'sampleSize':10})
    
    ####################   big    ##################################################
    
#     computeDatasetProbs({'modelName':'PubLDGenoGuard',
#                          'genotypeFileName':'./hapmap/new_genotypes_chr22_CEU.txt',
#                          'AFFileName':'./hapmap/hg_allele_freqs_chr22_CEU.txt',
#                          'LDFileName':'./hapmap/hg_ld_chr22_CEU.txt',
#                          'outFileName':'./tmp/publd_dataset0.txt',
#                          'idxBegin':0,
#                          'idxEnd':10})
#     generateModelSample({'modelName':'PubLDGenoGuard',
#                         'genotypeFileName':'./hapmap/new_genotypes_chr22_CEU.txt',
#                         'AFFileName':'./hapmap/hg_allele_freqs_chr22_CEU.txt',
#                         'LDFileName':'./hapmap/hg_ld_chr22_CEU.txt',
#                         'outFileName':'./tmp/publd_sample0.txt',
#                         'sampleSize':10})
#     computeDatasetProbs({'modelName':'DirectCondProbGenoGuard',
#                          'genotypeFileName':'./hapmap/new_genotypes_chr22_CEU.txt',
#                          'probFileName':'./hapmap/imm2_new_chr22_CEU_ref.txt',
#                          'outFileName':'./tmp/directcondprob_dataset2.txt',
#                          'order':2,
#                          'idxBegin':0,
#                          'idxEnd':10})
#     generateModelSample({'modelName':'DirectCondProbGenoGuard',
#                         'genotypeFileName':'./hapmap/new_genotypes_chr22_CEU.txt',
#                         'probFileName':'./hapmap/imm2_new_chr22_CEU_ref.txt',
#                         'outFileName':'./tmp/directcondprob_sample2.txt',
#                         'order':2,
#                         'sampleSize':10})
#     computeDatasetProbs({'modelName':'DirectCondProbGenoGuard',
#                          'genotypeFileName':'./hapmap/new_genotypes_chr22_CEU.txt',
#                          'probFileName':'./hapmap/imm1_new_chr22_CEU_ref.txt',
#                          'outFileName':'./tmp/directcondprob_dataset1.txt',
#                          'order':1,
#                          'idxBegin':0,
#                          'idxEnd':10})
#     generateModelSample({'modelName':'DirectCondProbGenoGuard',
#                         'genotypeFileName':'./hapmap/new_genotypes_chr22_CEU.txt',
#                         'probFileName':'./hapmap/imm1_new_chr22_CEU_ref.txt',
#                         'outFileName':'./tmp/directcondprob_sample1.txt',
#                         'order':1,
#                         'sampleSize':10})
#     computeDatasetProbs({'modelName':'DirectCondProbGenoGuard',
#                          'genotypeFileName':'./hapmap/new_genotypes_chr22_CEU.txt',
#                          'probFileName':'./hapmap/imm0_new_chr22_CEU_ref.txt',
#                          'outFileName':'./tmp/directcondprob_dataset0.txt',
#                          'order':0,
#                          'idxBegin':0,
#                          'idxEnd':10})
#     generateModelSample({'modelName':'DirectCondProbGenoGuard',
#                         'genotypeFileName':'./hapmap/new_genotypes_chr22_CEU.txt',
#                         'probFileName':'./hapmap/imm0_new_chr22_CEU_ref.txt',
#                         'outFileName':'./tmp/directcondprob_sample0.txt',
#                         'order':0,
#                         'sampleSize':10})