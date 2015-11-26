'''
Created on Mar 7, 2015

@author: zhihuang
'''
import random
import matplotlib.pyplot as plt
import numpy as np

def zipf(C, r, s):
    return C/(r**s)

def prob():
    total = 15250838
    unique = 486118
    C = 640505
    s = 0.905773
    p_ret = 0.699
    
    minRank = int(unique*(1-p_ret))
    f = sum([zipf(C, r, s) for r in range(minRank, unique)])
    total = sum([zipf(C, r, s) for r in range(1, unique+1)])
    P1 = zipf(C, 1, s)
    relative_P1 = P1*1.0 / (P1+f)
    ATT = 0.939
    worst_adv = ATT*relative_P1
    
    avg_adv = P1*1.0 / (P1 + (unique*p_ret-1)/unique*(total-P1))
    print "minRank: ", minRank, ", f: ", f, ", total: ", total, ", prob: ", f*1.0/total, "W: ", C*1.0/total, "relative_P1: ", relative_P1, "worst_adv: ", worst_adv, "avg_adv: ", avg_adv
    print 0.037871 / (0.037871+0.147651)
    print zipf(C, 2, s) / (total - P1)
    
def adv_analysis():
    total = 15250838
    unique = 486118
    C = 640505
    s = 0.905773
    p_ret = 0.699
    
    priors = [0.088, 0.426, 0.393, 0.093]
    accuracies = [[0.607, 0.286, 0.071, 0.036],
                [0.008, 0.939, 0.038, 0.015],
                [0.008, 0.567, 0.2, 0.225],
                [0, 0.552, 0.034, 0.414]]
    y1 = [zipf(C, r, s) for r in range(1, unique+1)]
    
    advs = [[]]
    rescaledMaxWeights = [[]]
    numExps = 1
    k = 0
    while(k < numExps):
        advs.append([])
        rescaledMaxWeights.append([])
        for i in range(4):
            retainedProb = 0
            for j in range(4):
                retainedProb += priors[j] * accuracies[j][i]
            print retainedProb
            y2 = [x for x in y1]
            for j in range(len(y2)):
                r = random.random()
                if  r > retainedProb:
                    y2[j] = 0
            rescaledMaxWeight = max(y2) / sum(y2)
            adv = accuracies[i][i] * rescaledMaxWeight
            advs[k].append(adv)
            rescaledMaxWeights[k].append(rescaledMaxWeight)
            print "accuracy: ", accuracies[i][i], ", rescaledMaxWeight: ", rescaledMaxWeight, ", adv: ", adv, ", originalMaxWeight: ", y1[0]/sum(y1)
        k += 1
    print "------------Average----------------"
    print "------------Red--------------------"
    avgRescaledMaxWeight = sum([rescaledMaxWeights[i][0] for i in range(numExps)]) / numExps
    avgAdv = sum([advs[i][0] for i in range(numExps)]) / numExps
    print "accuracy: ", accuracies[0][0], ", avgRescaledMaxWeight: ", avgRescaledMaxWeight, ", avgAdv: ", avgAdv, ", originalMaxWeight: ", y1[0]/sum(y1)
    print "------------Blond------------------"
    avgRescaledMaxWeight = sum([rescaledMaxWeights[i][1] for i in range(numExps)]) / numExps
    avgAdv = sum([advs[i][1] for i in range(numExps)]) / numExps
    print "accuracy: ", accuracies[1][1], ", avgRescaledMaxWeight: ", avgRescaledMaxWeight, ", avgAdv: ", avgAdv, ", originalMaxWeight: ", y1[0]/sum(y1)
    print "------------Brown------------------"
    avgRescaledMaxWeight = sum([rescaledMaxWeights[i][2] for i in range(numExps)]) / numExps
    avgAdv = sum([advs[i][2] for i in range(numExps)]) / numExps
    print "accuracy: ", accuracies[2][2], ", avgRescaledMaxWeight: ", avgRescaledMaxWeight, ", avgAdv: ", avgAdv, ", originalMaxWeight: ", y1[0]/sum(y1)
    print "------------Black------------------"
    avgRescaledMaxWeight = sum([rescaledMaxWeights[i][3] for i in range(numExps)]) / numExps
    avgAdv = sum([advs[i][3] for i in range(numExps)]) / numExps
    print "accuracy: ", accuracies[3][3], ", avgRescaledMaxWeight: ", avgRescaledMaxWeight, ", avgAdv: ", avgAdv, ", originalMaxWeight: ", y1[0]/sum(y1)
    
def plotAdv(orgAdv, newAdvs):
    fig, ax = plt.subplots()
    plt.rc('text', usetex=True)
    
    n = len(newAdvs)
    
    orgAdvs = [orgAdv] * n
    ind = np.arange(1, n+1)
    width = 0.35
    rects1 = ax.bar(ind, orgAdvs, width, color='r')
    rects2 = ax.bar(ind+width, newAdvs, width, color='g')
    
    ax.set_xlabel(r"Victim's Hair Color Known by Adversary $\mathcal{B}^\prime$", fontsize=14)
    ax.set_ylabel('Adversary\'s Advantage', fontsize=14)
    ax.set_xticks(ind+width)
    ax.set_xticklabels(('Red', 'Blond', 'Brown', 'Black'))
    
    ax.legend((rects1[0], rects2[0]), (r'Adversary $\mathcal{B}$', r'Adversary $\mathcal{B}^\prime$'))
    
    def autolabel(rects):
        for rect in rects:
            height = rect.get_height()
            ax.text(rect.get_x()+rect.get_width()/2., 1.03*height, '%.4f' % (height), ha = 'center', va='bottom')
    
    autolabel(rects1)
    autolabel(rects2)
    
    plt.show()
    
if __name__ == '__main__':
    prob()
#     adv_analysis()
#    plotAdv(0.0378712253105, [0.0641933884619, 0.0415538119276, 0.0188922866065, 0.038193388665])