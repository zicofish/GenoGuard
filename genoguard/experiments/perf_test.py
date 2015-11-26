#!/usr/bin/python
'''
Created on Nov 13, 2014

@author: zhihuang
'''
from honey_encryption.GenoGuard import *

def cluster_encode(kwargs):
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
    
    haplotype, geneticDist = loadRecombData(kwargs['haplotypeFileName'], kwargs['geneticDistFileName'])
    test_gg = RecombGenoGuard(SNPRefList, bitLen, map(lambda u: u[:10], haplotype), geneticDist)
    perf = test_gg.perfTest(dataset[:, int(kwargs['seqIdx'])], 'password')
    perfFile = open(kwargs['perfFileName'], 'w')
    perfFile.write(' '.join(map(str, perf)))
    perfFile.close()
    
def plotPerfResult():
    chromLens = []
    for i in range(1, 23):
        datFile = open('./hapmap/chr%s/big_genotypes_chr%s_CEU.txt' % (i, i))
        chromLens.append(len(datFile.readlines()))
        datFile.close()
    print max(chromLens), ' ', min(chromLens)
    
    chromPerf = []
    for i in range(1, 23):
        perfFile = open('./merge-output/hapmap/chr%s/big_recomb_perf_chr%s_CEU.txt' % (i, i))
        perfs = map(lambda u: map(float, u.split()), perfFile.readlines())
        perfs = np.array(perfs)
        avg_perf = np.average(perfs, axis = 0)
        chromPerf.append(avg_perf)
    chromPerf = np.array(chromPerf)
    
    import matplotlib.pyplot as plt
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    h1 = ax.plot(range(1, 23), list(chromPerf[:, 0]), marker='s', c='green', label='encode')
    h2 = ax.plot(range(1, 23), list(chromPerf[:, 1]), marker='^', c='red', label='decode')
    h3 = ax.plot(range(1, 23), list(chromPerf[:, 2]), marker='*', c='blue', label='PBE_encrypt')
    h4 = ax.plot(range(1, 23), list(chromPerf[:, 3]), marker='+', c='Magenta', label='PBE_decrypt')
    ax.set_xlabel('Chromosome', fontsize=16)
    ax.set_ylabel(r'Time (seconds)', fontsize=16)
    ax.legend(bbox_to_anchor=(0.998, 0.995),
              fontsize=12)
    plt.xlim(0, 23)
    plt.ylim(-1, 55)
    
    ax2 = ax.twinx()
    ax2.plot(range(1, 23), chromLens, marker='o', linestyle='--', c='black', label='Chromosome length')
    ax2.set_ylabel(r'Chromosome length (# of SNVs)', fontsize=16)
    ax2.legend(bbox_to_anchor=(0.998, 0.74),
              fontsize=12)
    
    plt.xticks(range(1, 23), map(lambda u: '%s' % (u), range(1, 23)))
    
    plt.show()
    
if __name__ == '__main__':
    #####  cluster only   ########
#     import sys
#     arg = sys.argv
#     arg_dict = dict(map(lambda u: u.split('='), arg[2:]))
#     current_module = sys.modules[__name__]
#     getattr(current_module, arg[1])(arg_dict)
    #####  cluster end   #########
    plotPerfResult()
#     cluster_encode({'modelName':'RecombGenoGuard',
#                      'genotypeFileName':'./hapmap/chr%s/big_genotypes_chr%s_CEU.txt' % (22, 22),
#                      'haplotypeFileName':'./hapmap/chr%s/big_CEU.chr%s.hap' % (22, 22),
#                      'geneticDistFileName':'./hapmap/chr%s/big_genetic_map_chr%s_combined_b36.txt' % (22, 22),
#                      'seqIdx':0,
#                      'perfFileName':'./hapmap/chr%s/big_recomb_perf_chr%s_CEU.txt'})