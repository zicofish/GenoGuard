'''
Created on Jul 16, 2014

@author: zhihuang
'''
from honey_encryption.hg import *
import gmpy2
from gmpy2 import mpz

def intervalSize(seq, SNPRefList, AF, LD, bitLen):
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
    return U[n%2] - L[n%2] + 1

def intervalSize2(seed, SNPRefList, AF, LD, bitLen = 1000):
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
    return U[n%2] - L[n%2] + 1

def naiveDecode(seed, SNPRefList, AF, LD, maxSeed):
    L = [-1, -1]; U = [-1, -1]
    L[0] = mpz(0); U[0] = mpz(maxSeed)
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

def plotResult(resultFileName):
    resultFile = open(resultFileName)
    correctValue = float(resultFile.readline())
    correctPwd = random.randint(1, 1001)
    naiveResult = map(float, resultFile.readline().split())
    honeyResult = map(float, resultFile.readline().split())
#     colors = ['blue'] * len(naiveResult)
#     colors[orgInd] = 'red'
#     sizes = [50] * len(naiveResult)
#     sizes[orgInd] = 150
    pwds = range(1, 1002)
    del pwds[correctPwd-1]
    
    import matplotlib.pyplot as plt
    fig = plt.figure()
    ax1 = fig.add_subplot(1, 2, 1)
    h1 = ax1.scatter(pwds, naiveResult, marker='o', s = 50, c='blue')
    h2 = ax1.scatter([correctPwd], [correctValue], marker='<', s = 500, c = 'red')
    ax1.set_xlabel('Password', fontsize=16)
    ax1.set_ylabel('Logarithm of Interval Size', fontsize=16)
    ax1.legend((h1, h2),
           ('Wrong Sequence', 'Correct Sequence'),
           bbox_to_anchor=(0.50, 1.02, 1., .102),
           scatterpoints=1,
           loc='lower right',
           ncol=2,
           fontsize=12)
    
    ax2 = fig.add_subplot(1, 2, 2)
    ax2.scatter(pwds, honeyResult, marker='o', s = 50, c='blue')
    ax2.scatter([correctPwd], [correctValue], marker='<', s = 500, c = 'red')
    ax2.set_xlabel('Password', fontsize=16)
    ax2.set_ylabel('Logarithm of Interval Size', fontsize=16)
    
#    m_corrects = [38032.885838797571, 38293.300203513041, 38678.328355820006, 37735.715451576172, 37787.591629935734, 38660.700367593316, 38607.948218205325, 37798.902131875759, 38075.408862743228, 38275.875517068227, 38087.55752857749, 37746.09741144177, 37573.108634785211, 38026.612754665432, 37647.277382773485, 38010.4351788899, 38351.943913870993, 38027.82979003682, 38924.858392380571, 37579.177483982348, 37914.814177640059, 37890.041684974138, 38409.915552435035, 38276.426622527979, 37984.26446384078, 38294.293941655669, 38453.4207699299, 38049.549772735729, 38350.415641990294, 38729.065320255046, 38141.99239923514, 38380.175461237377, 38237.361100540067, 38443.791582113627, 38241.889053953106, 37627.960192235274, 37535.059486288956, 37129.799719317351, 37919.073854344148, 37494.169004905147, 37324.722901071378, 38194.204751658413, 38100.361774353987, 38410.157654695569, 39147.176189643535, 38711.946134383026, 38253.152977758335, 37888.585978870338, 38434.331842913904, 37467.93596746424]
#    ax2.scatter([695]*50, m_corrects, marker='<', s = 500, c = 'red')
    
    plt.show()
    
def plotResultForOneVictim(resultFileName):
    resultFile = open(resultFileName)
    correctIdx = int(resultFile.readline())
    naiveResult = map(float, resultFile.readline().split())
    honeyResult = map(float, resultFile.readline().split())
#     colors = ['blue'] * len(naiveResult)
#     colors[orgInd] = 'red'
#     sizes = [50] * len(naiveResult)
#     sizes[orgInd] = 150
    pwds = range(1, 101)
    del pwds[correctIdx]
    correctValue = naiveResult[correctIdx]
    del naiveResult[correctIdx]
    del honeyResult[correctIdx]
    
    import matplotlib.pyplot as plt
    fig = plt.figure()
    ax1 = fig.add_subplot(1, 2, 1)
    h1 = ax1.scatter(pwds, naiveResult, marker='o', s = 50, c='blue')
    h2 = ax1.scatter(pwds[correctIdx:(correctIdx+1)], [correctValue], marker='<', s = 500, c = 'red')
    ax1.set_xlabel('Password', fontsize=16)
    ax1.set_ylabel('Logarithm of Interval Size', fontsize=16)
    ax1.legend((h1, h2),
           ('Wrong Sequence', 'Correct Sequence'),
           bbox_to_anchor=(0.50, 1.02, 1., .102),
           scatterpoints=1,
           loc='lower right',
           ncol=2,
           fontsize=12)
    
    ax2 = fig.add_subplot(1, 2, 2)
    ax2.scatter(pwds, honeyResult, marker='o', s = 50, c='blue')
    ax2.scatter(pwds[correctIdx:(correctIdx+1)], [correctValue], marker='<', s = 500, c = 'red')
    ax2.set_xlabel('Password', fontsize=16)
    ax2.set_ylabel('Logarithm of Interval Size', fontsize=16)
    
    plt.show()
    
    
def produce_results():
    test_AF, test_LD = loadPriorKnowledge("E:/Lab/LCA1/honeygenes/dat/hg_allele_freqs_chr22_CEU.txt",
                                "E:/Lab/LCA1/honeygenes/dat/hg_ld_chr22_CEU.txt")
    
    SNPRefList = []
    dataset = []
    datFile = open("E:/Lab/LCA1/honeygenes/dat/new_genotypes_chr22_CEU.txt")
    for line in datFile.readlines():
        attrArray = line.split()
        SNPRefList.append(attrArray[0])
        dataset.append(attrArray[1:])
#    SNPRefList = SNPRefList[0:100]
    dataset = np.array(dataset, dtype=int)
#    bitLen1 = int(math.ceil(len(SNPRefList)*math.log(3,2)))
    maxNaiveSeed = 3**len(SNPRefList) - 1
    bitLen2 = int(math.ceil(len(SNPRefList)*3))
    for i in range(len(SNPRefList)):
        powOfThree.append(mpz(3) ** mpz(i))

    
    for i in range(50):
        print gmpy2.log2(intervalSize(list(dataset[:, i]), SNPRefList, test_AF, test_LD, bitLen2))
#    sec_exp_result = open('sec_exp_roc_result_large3.txt', 'w')
#    pInd = 16
#    while pInd < 17:
#        pInd += 1
#        correctSize = gmpy2.log2(intervalSize(list(dataset[:, pInd]), SNPRefList, test_AF, test_LD, bitLen2))
#        numOfExps = 500
#        sizes1 = [0] * numOfExps
#        for i in range(numOfExps):
#            randomSeed1 = random.randint(0, maxNaiveSeed)
#            randomSeq1 = naiveDecode(randomSeed1, SNPRefList, test_AF, test_LD, maxNaiveSeed)
#            sizes1[i] = gmpy2.log2(intervalSize(randomSeq1, SNPRefList, test_AF, test_LD, bitLen2))
#        print ', '.join(map(str, sizes1))
#           
#        sizes2 = [0] * numOfExps
#        for i in range(numOfExps):
#            randomSeed2 = random.randint(0, 2**bitLen2 - 1)
#            sizes2[i] = gmpy2.log2(intervalSize2(randomSeed2, SNPRefList, test_AF, test_LD, bitLen2))
#        print ', '.join(map(str, sizes2))
#        sec_exp_result.write(str(correctSize)+'\n')
#        sec_exp_result.write('\t'.join(map(str, sizes1)) + '\n')
#        sec_exp_result.write('\t'.join(map(str, sizes2)) + '\n')
#    sec_exp_result.close()

def eval_models():
    test_AF, test_LD = loadPriorKnowledge("E:/Lab/LCA1/honeygenes/dat/hg_allele_freqs_chr22_CEU.txt",
                                "E:/Lab/LCA1/honeygenes/dat/hg_ld_chr22_CEU.txt")
    
    SNPRefList = []
    dataset = []
    datFile = open("E:/Lab/LCA1/honeygenes/dat/new_genotypes_chr22_CEU.txt")
    for line in datFile.readlines():
        attrArray = line.split()
        SNPRefList.append(attrArray[0])
        dataset.append(attrArray[1:])
     
if __name__ == "__main__":
    plotResultForOneVictim("./tmp/sec_exp_result.txt")
#    produce_results()
#    plotResult("sec_exp_roc_result_large30.txt")
#     test_AF, test_LD = loadPriorKnowledge("C:/Users/zhihuang/Desktop/dat/hapmap/hg_allele_freqs_chr22_CEU.txt",
#                                 "C:/Users/zhihuang/Desktop/dat/hapmap/hg_ld_chr22_CEU.txt")
#     
#     SNPRefList = []
#     dataset = []
#     datFile = open("C:/Users/zhihuang/Desktop/dat/hapmap/new_genotypes_chr22_CEU.txt")
#     for line in datFile.readlines():
#         attrArray = line.split()
#         SNPRefList.append(attrArray[0])
#         dataset.append(attrArray[1:])
# #     SNPRefList = SNPRefList[0:10000]
#     dataset = np.array(dataset, dtype=int)
#     bitLen1 = int(math.ceil(len(SNPRefList)*math.log(3,2)))
#     bitLen2 = int(math.ceil(len(SNPRefList)*4))
#     print bitLen1, bitLen2
#     for i in range(len(SNPRefList)):
#         powOfThree.append(mpz(3) ** mpz(i))
#     
#     numOfExps = 100
#     sizes1 = [0] * numOfExps
#     correct = random.randint(0, numOfExps - 1)
#     sizes1[correct] = gmpy2.log2(intervalSize(list(dataset[:, 17]), SNPRefList, test_AF, test_LD, bitLen2))
#     for i in range(numOfExps):
#         if i == correct: continue
#         randomSeed1 = random.randint(0, 2**bitLen1-1)
#         randomSeq1 = decode(randomSeed1, SNPRefList, test_AF, test_LD, bitLen1)
#         sizes1[i] = gmpy2.log2(intervalSize(randomSeq1, SNPRefList, test_AF, test_LD, bitLen2))
#     print correct
#     print ', '.join(map(str, sizes1))
#     
#     sizes2 = [0] * numOfExps
#     sizes2[correct] = sizes1[correct]
#     for i in range(numOfExps):
#         if i == correct: continue
#         randomSeed2 = random.randint(0, 2**bitLen2 - 1)
#         randomSeq2 = decode(randomSeed2, SNPRefList, test_AF, test_LD, bitLen2)
#         sizes2[i] = gmpy2.log2(intervalSize(randomSeq2, SNPRefList, test_AF, test_LD, bitLen2))
#     print ', '.join(map(str, sizes2))
#     
#     sec_exp_result = open('sec_exp_result.txt', 'w')
#     sec_exp_result.write(str(correct)+'\n')
#     sec_exp_result.write('\t'.join(map(str, sizes1)) + '\n')
#     sec_exp_result.write('\t'.join(map(str, sizes2)) + '\n')
#     sec_exp_result.close()
    
    
#     randomSeed2 = random.randint(0, 2**bitLen2 - 1)
#     randomSeq2 = decode(randomSeed2, SNPRefList, test_AF, test_LD, bitLen2)
#     size2 = intervalSize(randomSeq2, SNPRefList, test_AF, test_LD, bitLen2)
#     print gmpy2.log2(size2)
# #     
# #     ratio = size2 / size1
# #     print gmpy2.log2(ratio)
#     for i in range(np.shape(dataset)[1]):
#         tmp = gmpy2.log2(intervalSize(list(dataset[:, i]), SNPRefList, test_AF, test_LD, bitLen2))
#         print tmp