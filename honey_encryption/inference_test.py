'''
Created on Dec 22, 2014

@author: zhihuang
'''
def estimateError(kwargs):
    genotypeFile = open(kwargs['genotypeFileName'])
    genotype = map(lambda u: u.split()[int(kwargs['stardIdx']):int(kwargs['endIdx'])], genotypeFile.readlines())
    
    predictedGenotypeFile = open(kwargs['predictedGenotypeFileName'])
    predictedGenotype = map(lambda u: u.split(), predictedGenotypeFile.readlines())
    
    acc = 0
    for i in range(len(genotype)):
        for j in range(len(genotype[0])):
            acc += abs(int(genotype[i][j])-int(predictedGenotype[i][j]))
    return acc*1.0/(len(genotype)*len(genotype[0])*kwargs['hidPercent'])

def plotInferResult(hiddenPercentages, errors):
    import matplotlib.pyplot as plt
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    h1 = ax.plot(hiddenPercentages, errors, marker='o', c='green', label='encode')
    ax.set_xlabel('hidden percentage', fontsize=16)
    ax.set_ylabel(r'Average Estimation Error', fontsize=16)
#     ax.legend(bbox_to_anchor=(0.998, 0.995),
#               fontsize=12)
#     plt.xlim(0, 23)
#     plt.ylim(-1, 55)
    
#     ax2 = ax.twinx()
#     ax2.plot(range(1, 23), chromLens, marker='o', linestyle='--', c='black', label='Chromosome length')
#     ax2.set_ylabel(r'Chromosome length (# of SNVs)', fontsize=16)
#     ax2.legend(bbox_to_anchor=(0.998, 0.74),
#               fontsize=12)
#     
#     plt.xticks(range(1, 23), map(lambda u: '%s' % (u), range(1, 23)))
    
    plt.show()

if __name__ == '__main__':
    errors = []
    for i in range(1, 21):
        e = estimateError({'genotypeFileName':'./hapmap/chr22/big_genotypes_chr22_CEU.txt',
                   'stardIdx':100,
                   'endIdx':165,
                   'hidPercent':i*0.05,
                   'predictedGenotypeFileName':'./merge-output/hapmap/chr22/infer_test/PredictedSNVs_chr22_CEU_%s.txt' % (i*0.05)})
        errors.append(e)
    plotInferResult([i*0.05 for i in range(1, 21)], errors)