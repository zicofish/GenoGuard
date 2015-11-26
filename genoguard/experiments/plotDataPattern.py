'''
Created on Dec 11, 2014

@author: zhihuang
'''
from numpy import *

def loadDataSet(fileName):
    fr = open(fileName)
    datArr = map(lambda u: map(int, u.split())[1:], fr.readlines())
    return mat(datArr).transpose()

def pca(dataMat, topNfeat = 9999999):
    meanVals = mean(dataMat, axis = 0)
    meanRemoved = dataMat - meanVals
    covMat = cov(meanRemoved, rowvar = 0)
    eigVals, eigVects = linalg.eig(mat(covMat))
    eigValInd = argsort(eigVals)
    eigValInd = eigValInd[:-(topNfeat+1):-1]
    redEigVects = eigVects[:, eigValInd]
    lowDDataMat = meanRemoved * redEigVects
    reconMat = (lowDDataMat * redEigVects.T) + meanVals
    return redEigVects, lowDDataMat, reconMat

def plotPattern():
    dataMat = loadDataSet('./hapmap/chr22/big_genotypes_chr22_CEU.txt')
    eigVects, lowDMat, reconMat = pca(dataMat[:, 1:1000], 2)
    
    import matplotlib.pyplot as plt
    fig = plt.figure()
    
    ax1 = fig.add_subplot(2, 2, 1)

    ax1.scatter(lowDMat[:, 0].flatten().A[0], lowDMat[:, 1].flatten().A[0])
    
    dataMat = loadDataSet('./hapmap/chr22/pca_big_directcond_random_chr22_CEU.txt')
    eigVects, lowDMat, reconMat = pca(dataMat[:, 1:1000], 2)

    ax1.scatter(lowDMat[:, 0].flatten().A[0], lowDMat[:, 1].flatten().A[0], c='r')
    
    dataMat = loadDataSet('./hapmap/chr22/synthetic_recomb_chr22_CEU.txt')
    eigVects, lowDMat, reconMat = pca(dataMat[:, 1:1000], 2)

    ax1.scatter(lowDMat[:, 0].flatten().A[0], lowDMat[:, 1].flatten().A[0], c='g')
    plt.show()

if __name__ == "__main__":
    plotPattern()