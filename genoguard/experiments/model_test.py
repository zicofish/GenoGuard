#!/usr/bin/python
'''
Created on Nov 7, 2014

@author: zhihuang

Chi-squared test of different data models. Refer to the paper for more details.
'''
from honey_encryption.GenoGuard import *
import numpy as np
from scipy import stats
import pickle

def calcBins(maxSeed, k = 10):
    bins = []
    for i in range(k):
        bins.append(maxSeed*(i+1)/k)
    return bins

def calcModelBins(ggModel):
    bins = calcBins(ggModel.getMaxSeed())
    modelBins = []
    for i in range(len(bins)):
        seq = ggModel.decode(bins[i])
        modelBins.append(gmpy2.mpz(''.join(map(str, seq)), 3))
    return modelBins

def loadModelBins(binSeqFileName):
    binSeqFile = open(binSeqFileName)
    modelBins = map(lambda u: gmpy2.mpz(''.join(u.split()), 3), binSeqFile.readlines())
    return modelBins

def clusterCalcModelBins(kwargs):
    SNPRefList = []
    dataset = []
    datFile = open(kwargs['genotypeFileName'])
    for line in datFile.readlines():
        attrArray = line.split()
        SNPRefList.append(attrArray[0])
        dataset.append(attrArray[1:])
    datFile.close()
    dataset = np.array(dataset, dtype=int)
    bitLen = int(math.ceil(len(SNPRefList)*4))    # 3.5 is the storage overhead parameter as introduced in the paper.
    bitLen = (bitLen / 128 + 1) * 128    # Just for AES encryption convenience
    
    if kwargs['modelName'] == 'PubLDGenoGuard':
        test_AF, test_LD = loadPriorKnowledge(kwargs['AFFileName'], kwargs['LDFileName'])
        test_gg = PubLDGenoGuard(SNPRefList, bitLen, test_AF, test_LD)
    elif kwargs['modelName'] == 'DirectCondProbGenoGuard':
        directCondProbs = loadDirectCondProbs(kwargs['probFileName'])
        test_gg = DirectCondProbGenoGuard(SNPRefList, bitLen, directCondProbs, int(kwargs['order']))
    elif kwargs['modelName'] == 'RecombGenoGuard':
        haplotype, geneticDist = loadRecombData(kwargs['haplotypeFileName'], kwargs['geneticDistFileName'])
        test_gg = RecombGenoGuard(SNPRefList, bitLen, map(lambda u: u[:100], haplotype), geneticDist)
    elif kwargs['modelName'] == 'NaiveGenoGuard':
        test_gg = NaiveGenoGuard(SNPRefList)
        
    bins = calcBins(test_gg.getMaxSeed())
    seq = test_gg.decode(bins[int(kwargs['binIdx'])])
    binSeqFile = open(kwargs['binSeqFileName'], 'w')
    binSeqFile.write(' '.join(map(str, seq)) + '\n')
    binSeqFile.close()
    

def modelTest(kwargs):
#     SNPRefList = []
#     dataset = []
#     datFile = open(kwargs['genotypeFileName'])
#     for line in datFile.readlines():
#         attrArray = line.split()
#         SNPRefList.append(attrArray[0])
#         dataset.append(attrArray[1:])
#     datFile.close()
#     dataset = np.array(dataset, dtype=int)
#     bitLen = int(math.ceil(len(SNPRefList)*4))    # 3.5 is the storage overhead parameter as introduced in the paper.
#     bitLen = (bitLen / 128 + 1) * 128    # Just for AES encryption convenience
#     
#     if kwargs['modelName'] == 'PubLDGenoGuard':
#         test_AF, test_LD = loadPriorKnowledge(kwargs['AFFileName'], kwargs['LDFileName'])
#         test_gg = PubLDGenoGuard(SNPRefList, bitLen, test_AF, test_LD)
#     elif kwargs['modelName'] == 'DirectCondProbGenoGuard':
#         directCondProbs = loadDirectCondProbs(kwargs['probFileName'])
#         test_gg = DirectCondProbGenoGuard(SNPRefList, bitLen, directCondProbs, int(kwargs['order']))
#     elif kwargs['modelName'] == 'RecombGenoGuard':
#         haplotype, geneticDist = loadRecombData(kwargs['haplotypeFileName'], kwargs['geneticDistFileName'])
#         test_gg = RecombGenoGuard(SNPRefList, bitLen, map(lambda u: u[:100], haplotype), geneticDist)
#     elif kwargs['modelName'] == 'NaiveGenoGuard':
#         test_gg = NaiveGenoGuard(SNPRefList)
    
#     modelBins = calcModelBins(test_gg)
    modelBins = loadModelBins(kwargs['binSeqFileName'])
    
    intSeqsFile = open(kwargs['intSeqsFileName'])
    intSeqs = pickle.load(intSeqsFile)
    intSeqs = sorted(intSeqs)
    intSeqsFile.close()
    
    observedFreq = []
    j = 0
    for i in range(len(modelBins)):
        observedFreq.append(0)
        while j < len(intSeqs):
            if intSeqs[j] < modelBins[i]:
                observedFreq[i] += 1
            else:
                break
            j += 1
    print observedFreq
    print stats.chisquare(observedFreq)
    return stats.chisquare(observedFreq)
    
def testAllChroms():
#     naiveStatFile = open('./hapmap/big_stat_naive.txt', 'w')
#     for i in range(1, 23):
#         stat = modelTest({'modelName':'NaiveGenoGuard',
#                         'genotypeFileName':'./hapmap/chr%s/big_genotypes_chr%s_CEU.txt' % (i, i),
#                         'intSeqsFileName':'./hapmap/chr%s/big_int_seqs_chr%s_CEU.pickle' % (i, i)})
#         naiveStatFile.write('\t'.join(map(str, stat)) + '\n')
#     naiveStatFile.close()
#     pubLDStatFile = open('./hapmap/big_stat_publd.txt', 'w')
#     for i in range(1, 23):
#         stat = modelTest({'modelName':'PubLDGenoGuard',
#                         'genotypeFileName':'./hapmap/chr%s/big_genotypes_chr%s_CEU.txt' % (i, i),
#                         'AFFileName':'./hapmap/chr%s/big_allele_freqs_chr%s_CEU.txt' % (i, i),
#                         'LDFileName':'./hapmap/chr%s/big_ld_chr%s_CEU.txt' % (i, i),
#                         'intSeqsFileName':'./hapmap/chr%s/big_int_seqs_chr%s_CEU.pickle' % (i, i)})
#         pubLDStatFile.write('\t'.join(map(str, stat)) + '\n')
#     pubLDStatFile.close()
#     condProb0StatFile = open('./hapmap/big_stat_condProb0.txt', 'w')
#     for i in range(1, 23):
#         stat = modelTest({'modelName':'DirectCondProbGenoGuard',
#                         'genotypeFileName':'./hapmap/chr%s/big_genotypes_chr%s_CEU.txt' % (i, i),
#                         'probFileName':'./hapmap/chr%s/big_condProb0_chr%s_CEU_ref.txt' % (i, i),
#                         'intSeqsFileName':'./hapmap/chr%s/big_int_seqs_chr%s_CEU.pickle' % (i, i),
#                         'order':0,
#                         'binSeqFileName': './merge-output/hapmap/chr%s/big_condProb0_chr%s_bin_seq_CEU.txt' % (i, i)})
#         condProb0StatFile.write('\t'.join(map(str, stat)) + '\n')
#   
#     condProb1StatFile = open('./hapmap/big_stat_condProb1.txt', 'w')
#     for i in range(1, 23):
#         stat = modelTest({'modelName':'DirectCondProbGenoGuard',
#                         'genotypeFileName':'./hapmap/chr%s/big_genotypes_chr%s_CEU.txt' % (i, i),
#                         'probFileName':'./hapmap/chr%s/big_condProb1_chr%s_CEU_ref.txt' % (i, i),
#                         'intSeqsFileName':'./hapmap/chr%s/big_int_seqs_chr%s_CEU.pickle' % (i, i),
#                         'order':1,
#                         'binSeqFileName': './merge-output/hapmap/chr%s/big_condProb1_chr%s_bin_seq_CEU.txt' % (i, i)})
#         condProb1StatFile.write('\t'.join(map(str, stat)) + '\n')
#     condProb1StatFile.close()
#   
#     condProb2StatFile = open('./hapmap/big_stat_condProb2.txt', 'w')
#     for i in range(1, 23):
#         stat = modelTest({'modelName':'DirectCondProbGenoGuard',
#                         'genotypeFileName':'./hapmap/chr%s/big_genotypes_chr%s_CEU.txt' % (i, i),
#                         'probFileName':'./hapmap/chr%s/big_condProb2_chr%s_CEU_ref.txt' % (i, i),
#                         'intSeqsFileName':'./hapmap/chr%s/big_int_seqs_chr%s_CEU.pickle' % (i, i),
#                         'order':2,
#                         'binSeqFileName': './merge-output/hapmap/chr%s/big_condProb2_chr%s_bin_seq_CEU.txt' % (i, i)})
#         condProb2StatFile.write('\t'.join(map(str, stat)) + '\n')
#     condProb2StatFile.close()
    recombStatFile = open('./hapmap/big_stat_recomb.txt', 'w')
    for i in range(1, 23):
        stat = modelTest({'modelName':'RecombGenoGuard',
                         'genotypeFileName':'./hapmap/chr%s/big_genotypes_chr%s_CEU.txt' % (i, i),
                         'haplotypeFileName':'./hapmap/chr%s/big_CEU.chr%s.hap' % (i, i),
                         'geneticDistFileName':'./hapmap/chr%s/big_genetic_map_chr%s_combined_b36.txt' % (i, i),
                         'intSeqsFileName':'./hapmap/chr%s/big_int_seqs_chr%s_CEU.pickle' % (i, i),
                         'binSeqFileName': './merge-output/hapmap/chr%s/big_recomb_chr%s_bin_seq_CEU.txt' % (i, i)})
        recombStatFile.write('\t'.join(map(str, stat)) + '\n')
    recombStatFile.close()
    
def plotTestResult():
    naiveStatFile = open('./hapmap/chi-square/big_stat_naive.txt')
    naiveStat = map(lambda u: map(float, u.split()), naiveStatFile.readlines())
    naiveStatFile.close()
    
    pubLDStatFile = open('./hapmap/chi-square/big_stat_publd.txt')
    pubLDStat = map(lambda u: map(float, u.split()), pubLDStatFile.readlines())
    pubLDStatFile.close()
    
    dataset0StatFile = open('./hapmap/chi-square/big_stat_condProb0.txt')
    dataset0Stat = map(lambda u: map(float, u.split()), dataset0StatFile.readlines())
    dataset0StatFile.close()
    
    dataset1StatFile = open('./hapmap/chi-square/big_stat_condProb1.txt')
    dataset1Stat = map(lambda u: map(float, u.split()), dataset1StatFile.readlines())
    dataset1StatFile.close()
    
    dataset2StatFile = open('./hapmap/chi-square/big_stat_condProb2.txt')
    dataset2Stat = map(lambda u: map(float, u.split()), dataset2StatFile.readlines())
    dataset2StatFile.close()
    
    recombStatFile = open('./hapmap/chi-square/big_stat_recomb.txt')
    recombStat = map(lambda u: map(float, u.split()), recombStatFile.readlines())
    recombStatFile.close()
    
    chi_thresholds = map(lambda u: math.log(stats.chi2.ppf(u, df=9)), [1-0.2, 1-0.01])
    
    import matplotlib.pyplot as plt
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    h1 = ax.plot(range(1, 23), map(lambda u: math.log(u[0]), naiveStat), marker='s', c='grey', label='Uniform distribution model')
    h2 = ax.plot(range(1, 23), map(lambda u: math.log(u[0]), pubLDStat), marker='o', c='red', label='Public LD model')
    h3 = ax.plot(range(1, 23), map(lambda u: math.log(u[0]), dataset0Stat), marker='*', c='blue', label='0-th-order model')
    h4 = ax.plot(range(1, 23), map(lambda u: math.log(u[0]), dataset1Stat), marker='+', c='Magenta', label='1-st-order model')
    h5 = ax.plot(range(1, 23), map(lambda u: math.log(u[0]), dataset2Stat), marker='^', c='Cyan', label='2-nd-order model')
    h6 = ax.plot(range(1, 23), map(lambda u: math.log(u[0]), recombStat), marker='D', c='green', label='Recombination model')
    h7 = ax.plot([0, 23], [chi_thresholds[0], chi_thresholds[0]], linestyle='--', c='black')
    h8 = ax.plot([0, 23], [chi_thresholds[1], chi_thresholds[1]], linestyle='--', c='black')
#     ax.text(23.3, math.log(chi_threshold)-0.1, r'$\alpha = 0.01$', fontsize = 18, weight='bold')
    ax.set_xlabel('Chromosome', fontsize=16)
    ax.set_ylabel(r'$\log(\chi^2)$', fontsize=16)
    ax.legend(loc='upper center',
              ncol=2,
              fontsize=12)
    plt.xlim(0, 23)
    plt.ylim(0, 9)
    plt.xticks(range(1, 23), map(lambda u: '%s' % (u), range(1, 23)))
    
    ax2 = ax.twinx()
    ax2.set_ylabel(r'Significance level $\alpha$', fontsize=16)
    plt.yticks(range(10) + chi_thresholds, ['']*10 + map(str, [0.2, 0.01]))
    
#     fig = plt.figure()
#     ax = fig.add_subplot(1,1,1)
#     h1 = ax.plot(range(1, 23), map(lambda u: -math.log10(u[1]), naiveStat), marker='s', c='grey', label='Conventional model')
#     h2 = ax.plot(range(1, 23), map(lambda u: -math.log10(u[1]), pubLDStat), marker='o', c='red', label='Public LD model')
#     h3 = ax.plot(range(1, 23), map(lambda u: -math.log10(u[1]), dataset0Stat), marker='*', c='blue', label='0-th-order model')
#     h4 = ax.plot(range(1, 23), map(lambda u: -math.log10(u[1]), dataset1Stat), marker='+', c='Magenta', label='1-st-order model')
#     h5 = ax.plot(range(1, 23), map(lambda u: -math.log10(u[1]), dataset2Stat), marker='^', c='Cyan', label='2-nd-order model')
#     h6 = ax.plot(range(1, 23), map(lambda u: -math.log10(u[1]), recombStat), marker='D', c='green', label='Recombination model')
#     h7 = ax.plot([0, 23], [-math.log10(0.01), -math.log10(0.01)], linestyle='--', c='black')
#     ax.text(23.3, -math.log10(0.01), r'$\alpha = 0.01$', fontsize = 18, weight='bold')
#     ax.set_xlabel('Chromosome', fontsize=16)
#     ax.set_ylabel(r'$-\log_{10}(p)$', fontsize=16)
#     ax.legend(loc='upper center',
#               ncol=2,
#               fontsize=12)
#     plt.xlim(0, 23)
#     plt.ylim(-15, 350)
#     plt.xticks(range(1, 23), map(lambda u: '%s' % (u), range(1, 23)))
#     plt.yticks(range(0, 350, 50), map(str, range(0, 350, 50)))
    plt.show()
    
if __name__ == '__main__':
    #####  cluster only   ########
#     import sys
#     arg = sys.argv
#     arg_dict = dict(map(lambda u: u.split('='), arg[2:]))
#     current_module = sys.modules[__name__]
#     getattr(current_module, arg[1])(arg_dict)
    #####  cluster end   #########
    
    plotTestResult()
#     testAllChroms()
#     clusterCalcModelBins(0, {'modelName':'NaiveGenoGuard',
#                              'genotypeFileName':'./hapmap/chr%s/big_genotypes_chr%s_CEU.txt' % (22, 22),
#                              'intSeqsFileName':'./hapmap/chr%s/big_int_seqs_chr%s_CEU.pickle' % (22, 22),
#                              'binSeqFileName':'./hapmap/chr%s/big_chr%s_bin_seq_%s_CEU.txt' % (22, 22, 0)})
#     (7.6666666666666679, 0.5680553672950962)
#     modelTest({'modelName':'PubLDGenoGuard',
#              'genotypeFileName':'./hapmap/small_genotypes_chr22_CEU.txt',
#              'AFFileName':'./hapmap/small_allele_freqs_chr22_CEU.txt',
#              'LDFileName':'./hapmap/small_ld_chr22_CEU.txt'})
#     modelTest({'modelName':'DirectCondProbGenoGuard',
#             'genotypeFileName':'./hapmap/small_genotypes_chr22_CEU.txt',
#             'probFileName':'./hapmap/condProb0_chr22_CEU_ref.txt',
#             'order':0})
#     modelTest({'modelName':'DirectCondProbGenoGuard',
#              'genotypeFileName':'./hapmap/small_genotypes_chr22_CEU.txt',
#              'probFileName':'./hapmap/condProb1_chr22_CEU_ref.txt',
#              'order':1})
#     modelTest({'modelName':'DirectCondProbGenoGuard',
#              'genotypeFileName':'./hapmap/small_genotypes_chr22_CEU.txt',
#              'probFileName':'./hapmap/condProb2_chr22_CEU_ref.txt',
#              'order':2})
#     modelTest({'modelName':'RecombGenoGuard',
#              'genotypeFileName':'./hapmap/small_genotypes_chr22_CEU.txt',
#              'haplotypeFileName':'./hapmap/haplotype/new_CEU.chr22.hap',
#              'geneticDistFileName':'./hapmap/haplotype/new_genetic_map_chr22_combined_b36.txt'})  