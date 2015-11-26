'''
Created on Jun 11, 2014

@author: zhihuang
'''
from numpy import *
def loadDataSet(fileName, delim='\t'):
    fr = open(fileName)
    stringArr = [line.strip().split(delim) for line in fr.readlines()]
    ancestryArr = [line[1] for line in stringArr if line[1] in ['CEU', 'ASW', 'CHB', 'TESTASW', 'TESTCEU', 'TESTCHB']]
    datArr = [map(float, line[2:]) for line in stringArr if line[1] in ['CEU', 'ASW', 'CHB', 'TESTASW', 'TESTCEU', 'TESTCHB']]
    return array(ancestryArr), mat(datArr)

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
    
def replaceNanWithMean():
    datMat = loadDataSet('secom.dat', ' ')
    numFeat = shape(datMat)[1]
    for i in range(numFeat):
        meanVal = mean(datMat[nonzero(~isnan(datMat[:, i].A))[0], i])
        datMat[nonzero(isnan(datMat[:, i].A))[0], i] = meanVal
    return datMat

def plotPCNum(eigVals):
    totalVariance = sum(eigVals)
    percentage = eigVals / totalVariance
    cumulative = cumsum(percentage)
    import matplotlib.pyplot as plt
    fig = plt.figure()
    ax1 = fig.add_subplot(121, xlabel="Principal Component ID", ylabel="Percentage of Variance")
    ax1.plot(array(range(1, 21)), percentage[:20], marker='^')
    
    ax2 = fig.add_subplot(122, xlabel="Principal Component Number", ylabel="Cumulative Variance")
    ax2.plot(array(range(1, 21)), cumulative[:20], marker = '^')
    
    plt.show()
    
def ancestryExp(ancestryList):
    ancestryArr, dataMat = loadDataSet('trainingSet.txt')
    eigVects, lowDMat, reconMat = pca(dataMat, 2)
    
    import matplotlib.pyplot as plt
    fig = plt.figure()
    subPlotTitles = {'Org':'(a)', 'ASW':'(b)', 'CEU':'(c)', 'CHB':'(d)'}
    
    ax1 = fig.add_subplot(2, 2, 1)
    ax1.set_xlim(-4, 8)
    ax1.text(1.6, -6.5, subPlotTitles['Org'])
#     ax1.scatter(lowDMat[:, 0].flatten().A[0], lowDMat[:, 1].flatten().A[0], 
#                marker='o', s = 50, 
#                c = [ancestryColor[ancestry] for ancestry in ancestryArr])
    h1 = scatterGroup(ax1, lowDMat, ancestryArr, 'ASW')
    h2 = scatterGroup(ax1, lowDMat, ancestryArr, 'CEU')
    h3 = scatterGroup(ax1, lowDMat, ancestryArr, 'CHB')
    ax1.legend((h1, h2, h3),
           ('ASW', 'CEU', 'CHB'),
           scatterpoints=1,
           loc='lower right',
           ncol=2,
           fontsize=8)
    
    for ancestry in ancestryList:
        test_ancestryArr, test_dataMat = loadDataSet('exp'+ancestry+'.txt')
        meanVals = mean(dataMat, axis = 0)
        test_meanRemoved = test_dataMat - meanVals
        test_lowDMat = test_meanRemoved * eigVects

        ax2 = fig.add_subplot(2, 2, ancestryList.index(ancestry)+2)
        ax2.set_xlim(-4, 8)
        ax2.text(1.6, -6.5, subPlotTitles[ancestry])
        h1 = scatterGroup(ax2, lowDMat, ancestryArr, 'ASW')
        h2 = scatterGroup(ax2, lowDMat, ancestryArr, 'CEU')
        h3 = scatterGroup(ax2, lowDMat, ancestryArr, 'CHB')
#         ax2.scatter(lowDMat[:, 0].flatten().A[0], lowDMat[:, 1].flatten().A[0], 
#                    marker='o', s = 50, 
#                    c = [ancestryColor[ancestry] for ancestry in ancestryArr])

#         ax2.scatter(test_lowDMat[:, 0].flatten().A[0], test_lowDMat[:, 1].flatten().A[0], 
#                    marker='+', s = 50, 
#                    c = ancestryColor['TEST'+ancestry], linewidths=1.5)
        h4 = scatterGroup(ax2, test_lowDMat, test_ancestryArr, 'TEST'+ancestry)
        ax2.legend((h1, h2, h3, h4),
           ('ASW', 'CEU', 'CHB', 'TEST_'+ancestry),
           scatterpoints=1,
           loc='lower right',
           ncol=2,
           fontsize=8)
        
    plt.show()
    
def scatterGroup(m_ax, dataMat, ancestryArr, ancestry):
    ancestryColor = {'ASW':'green', 'CEU':'blue', 'CHD':'lavender', 'GIH':'yellow', 
                     'CHB':'orange', 'JPT':'pink', 'LWK':'purple', 'MKK':'gray', 
                     'MEX':'cyan', 'TSI':'black', 'YRI':'magenta', 'TESTASW':'red', 'TESTCEU':'red', 'TESTCHB':'red'}
    ancestryMarker = {'ASW':'o', 'CEU':'o', 'CHB':'o', 'TESTASW':'+', 'TESTCEU':'+', 'TESTCHB':'+'}
    ancestryLinewidth = {'ASW':1, 'CEU':1, 'CHB':1, 'TESTASW':1.5, 'TESTCEU':1.5, 'TESTCHB':1.5}
    group = dataMat[nonzero(ancestryArr == ancestry)[0]]
    return m_ax.scatter(group[:, 0].flatten().A[0], group[:, 1].flatten().A[0], marker=ancestryMarker[ancestry], s = 50, c = ancestryColor[ancestry], linewidths=ancestryLinewidth[ancestry])
    
ancestryExp(['ASW', 'CEU', 'CHB'])