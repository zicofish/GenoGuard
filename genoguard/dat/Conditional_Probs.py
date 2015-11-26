'''
Created on 4 nov. 2014

@author: shariati
'''
import numpy as np
import numpy.matlib

# loading the training dataset
def loadTraningDataset(TrainingDatasetName):
    
    print ("LOADING DATASET")
    
    TrainingData = open(TrainingDatasetName)
    SNVs = []
    for line in TrainingData.readlines():
        attrArray = line.split()
        SNVs.append(attrArray[1:])
    TrainingData.close()
    
    return SNVs

def condNum_k_order(SNVs,k):
    
    condNum_k_order = []
    
    if k == 0:
        for x in SNVs:
            num_0 = float(x.count('0'))
            num_1 = float(x.count('1'))
            num_2 = float(x.count('2'))
            condNum_k_order.append([num_0,num_1,num_2])
            
    if k == 1:
        m = []
        for i in range (0,k+1):
            m.append([0,1,2])
            condNum_k_order = [[] for x in xrange(0,len(SNVs))]
        
        for i in range (0,k):
            for x in m[0]:
                for y in m[1]:
                    condNum_k_order[i].append(0)
        
        for i in range (k,len(SNVs)):
            for x in m[0]:
                for y in m[1]:
                    num = 0
                    for j in range (0,len(SNVs[0])):
                        if SNVs[i-1][j] == ('%s' %x) and SNVs[i][j] == ('%s' %y):
                            num = num+1
                    condNum_k_order[i].append(float(num))
    
    if k == 2:
        m = []
        for i in range (0,k+1):
            m.append([0,1,2])
            condNum_k_order = [[] for x in xrange(0,len(SNVs))]
        
        for i in range (0,k):
            for x in m[0]:
                for y in m[1]:
                    for t in m[2]:
                        condNum_k_order[i].append(0)
                        
        for i in range (k,len(SNVs)):
            for x in m[0]:
                for y in m[1]:
                    for t in m[2]:
                        num = 0
                        for j in range(0,len(SNVs[0])):
                            if SNVs[i-2][j] == ('%s' %x) and SNVs[i-1][j] == ('%s' %y) and SNVs[i][j] == ('%s' %t):
                                num = num+1
                        condNum_k_order[i].append(float(num))

    return condNum_k_order

# calculating the probability of SNV(i)
def condProb_0_order(SNVs):
    
    print ("CALCULATING 0-ORDER CONDITIONAL PROBABILITY")
    
    condProb_0 = []
    condNum_0 = condNum_k_order(SNVs, 0)
    for x in condNum_0:
        y = map(lambda x: x/len(SNVs[0]), x)
        condProb_0.append(y)
    
    return condProb_0

# calculating the probability of SNV(i) given SNV(i-1)
def condProb_1_order(SNVs, condProb_0):
    
    print ("CALCULATING 1-ORDER CONDITIONAL PROBABILITY")    
    
    condProb_1 = []
#     condProb_0 = condProb_0_order(SNVs)
    expand_condProb_0 = numpy.matlib.repmat(condProb_0, 1, 3)


    condNum_0 = condNum_k_order(SNVs, 0)
    expand_condNum_0 = np.repeat(condNum_0, 3, axis=1)
    expand_condNum_0[np.nonzero(expand_condNum_0 == 0)] = 1e-10
    
    condNum_1 = np.array(condNum_k_order(SNVs, 1))
    
    condProb_1.append(expand_condProb_0[0])
    
    tmp = condNum_1[1:len(condNum_1), :] / expand_condNum_0[0:len(expand_condNum_0)-1, :]
    
    condProb_1 = np.insert(tmp, 0, condProb_1, axis=0)
    
    return condProb_1  

# calculating the probability of SNV(i) given SNV(i-1) and SNV(i-2)
def condProb_2_order(SNVs, condProb_1):
    
    print ("CALCULATING 2-ORDER CONDITIONAL PROBABILITY")
    
    condProb_2 = []
#     condProb_1 = condProb_1_order(SNVs)
    expand_condProb_1 = numpy.matlib.repmat(condProb_1, 1, 3)
    
    condNum_1 = condNum_k_order(SNVs, 1)
    expand_condNum_1 = np.repeat(condNum_1, 3, axis=1)
    expand_condNum_1[np.nonzero(expand_condNum_1 == 0)] = 1e-10
    
    condNum_2 = np.array(condNum_k_order(SNVs, 2))
    
    condProb_2.append(expand_condProb_1[0])
    condProb_2.append(expand_condProb_1[1])

    
    tmp = condNum_2[1:len(condNum_2), :] / expand_condNum_1[0:len(expand_condNum_1)-1, :]
    
    condProb_2 = np.insert(tmp[1:len(tmp),:], 0, condProb_2, axis=0)
    
    return condProb_2

def buildAllChroms():
    for i in range(1, 23):
        # generating outputs
        buildChrom(i)
        
def buildChrom(ChromNum):
    snvs_in = loadTraningDataset("../hapmap/chr%s/small_genotypes_chr%s_CEU.txt" % (ChromNum, ChromNum))
        
    condProb0 = condProb_0_order(snvs_in)
    condProb1 = condProb_1_order(snvs_in, condProb0)
    condProb2 = condProb_2_order(snvs_in, condProb1)
    print ".................."
    condProb0File = open('../hapmap/chr%s/small_condProb0_chr%s_CEU_ref.txt' % (ChromNum, ChromNum), 'w')
    condProb0File.write("\n".join(map(lambda u: '\t'.join(map(str, u)), condProb0)))
    condProb0File.close()
    # print "\n".join(map(str, condProb0))
    print ".................."
    condProb1File = open('../hapmap/chr%s/small_condProb1_chr%s_CEU_ref.txt' % (ChromNum, ChromNum), 'w')
    condProb1File.write("\n".join(map(lambda u: '\t'.join(map(str, u)), condProb1)))
    condProb1File.close()
    # print "\n".join(map(str, condProb1))
    print ".................."
    condProb2File = open('../hapmap/chr%s/small_condProb2_chr%s_CEU_ref.txt' % (ChromNum, ChromNum), 'w')
    condProb2File.write("\n".join(map(lambda u: '\t'.join(map(str, u)), condProb2)))
    condProb2File.close()
    # print "\n".join(map(str, condProb2))
        
if __name__ == '__main__':
    #buildAllChroms()
    buildChrom(22)