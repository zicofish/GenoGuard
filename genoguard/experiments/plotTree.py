'''
Created on Jun 25, 2014

@author: zhihuang
'''
import matplotlib.pyplot as plt
import numpy as np
import random

decisionNode = dict(boxstyle = "square", fc = "0.8")
leafNode = dict(boxstyle="round", fc="0.8")
arrow_args = dict(arrowstyle="<|-")

def plotNode(nodeTxt, centerPt, parentPt, nodeType, myColor, myLinestyle):
    nodeType['color'] = myColor
    nodeType['linestyle'] = myLinestyle
    arrow_args['color'] = myColor
    arrow_args['linestyle'] = myLinestyle
    parentPt = (parentPt[0], parentPt[1] - 0.02)
    createPlot.ax1.annotate(nodeTxt, xy=parentPt, xycoords='axes fraction',
                            xytext=centerPt, textcoords='axes fraction',
                            va="center", ha="center", bbox=nodeType,
                            arrowprops=arrow_args, weight='bold', color = myColor)

def plotMidText(childPt, parentPt, txtString, myColor):
    xMid = (parentPt[0] - childPt[0])/2.0 + childPt[0]
    yMid = (parentPt[1] - childPt[1])/2.0 + childPt[1]
    createPlot.ax1.text(xMid, yMid, txtString, weight='bold', color = myColor)
    
def plotTree(myTree, parentPt, nodeTxt, path):
    numLeafs = getNumLeafs(myTree)
#     getTreeDepth(myTree)
    firstStr = myTree.keys()[0]
    childPt = (plotTree.xOff + (1.0 + float(numLeafs))/2.0/plotTree.totalW, plotTree.yOff)
    myColor = 'black'
    myLinestyle = 'solid'
#    if path=='' or path=='0' or path=='02' or path=='021' or path=='00' or path=='001' or path=='01' or path=='011':
#        myColor = 'red'
#        myLinestyle = 'dashed'
    if path=='' or path=='0' or path=='02' or path=='021':
        myColor = 'red'
        myLinestyle = 'dashed'
    plotMidText(childPt, parentPt, nodeTxt, myColor)
    plotNode(firstStr, childPt, parentPt, decisionNode, myColor, myLinestyle)
    secondDict = myTree[firstStr]
    plotTree.yOff = plotTree.yOff - 1.0/plotTree.totalD
    incr = 0
    for key in secondDict.keys():
        if type(secondDict[key]).__name__ == "dict":
            plotTree(secondDict[key], childPt, str(key), path+str(key))
        else:
            myColor = 'black'
            myLinestyle = 'solid'
#            if path+str(key) == '021' or path+str(key) == '001' or path+str(key) == '011':
#                myColor = 'red'
#                myLinestyle = 'dashed'
            if path+str(key) == '021':
                myColor = 'red'
                myLinestyle = 'dashed'
            plotTree.xOff = plotTree.xOff + 1.0/plotTree.totalW
            if path == '':
                plotNode(secondDict[key], (plotTree.xOff, plotTree.yOff), childPt, leafNode, myColor, myLinestyle)
            else:
                plotNode(secondDict[key], (plotTree.xOff, plotTree.yOff - incr), childPt, leafNode, myColor, myLinestyle)
            plotMidText((plotTree.xOff, plotTree.yOff), childPt, str(key), myColor)
        incr += 0.1
    plotTree.yOff = plotTree.yOff + 1.0/plotTree.totalD
    
def createPlot(inTree):
    fig = plt.figure(1, facecolor='white')
    fig.clf()
    createPlot.ax1 = plt.subplot(111, frameon=False)
    createPlot.ax1.axis('off')
    plotTree.totalW = float(getNumLeafs(inTree))
    plotTree.totalD = float(getTreeDepth(inTree)+1)
    plotTree.xOff = -0.5/plotTree.totalW
    plotTree.yOff = 1.0
    plotTree(inTree, (0.5, 1.0), '', '')
    
    plt.setp(createPlot.ax1.get_xticklabels(), visible=False)
    plt.setp(createPlot.ax1.get_yticklabels(), visible=False)
    plt.show()
    
def getNumLeafs(myTree):
    numLeafs = 0
    firstStr = myTree.keys()[0]
    secondDict = myTree[firstStr]
    for key in secondDict.keys():
        if type(secondDict[key]).__name__ == 'dict':
            numLeafs += getNumLeafs(secondDict[key])
        else:
            numLeafs += 1
    return numLeafs

def getTreeDepth(myTree):
    maxDepth = 0
    firstStr = myTree.keys()[0]
    secondDict = myTree[firstStr]
    for key in secondDict.keys():
        if type(secondDict[key]).__name__ == 'dict':
            thisDepth = 1 + getTreeDepth(secondDict[key])
        else:
            thisDepth = 1
        if thisDepth > maxDepth:
            maxDepth = thisDepth
    return maxDepth
    
randCond = [0.5, 0.3, 0.2, 0.6, 0.3, 0.1]
def createTree(interval, depth):
    key = ', '.join(map(str, interval))
    if depth == 3:
        return key
    myTree = {key : {}}
    pick = random.randint(0, 1)
    p0 = randCond[pick*3]; p1 = randCond[pick*3 + 1]; p2 = randCond[pick*3+2]
    myTree[key][0] = createTree([interval[0], interval[0] + (interval[1]-interval[0])*p0], depth+1)
    myTree[key][1] = createTree([interval[0] + (interval[1]-interval[0])*p0, interval[0] + (interval[1] - interval[0]) * (p0+p1)], depth+1)
    myTree[key][2] = createTree([interval[0] + (interval[1] - interval[0]) * (p0+p1), interval[1]], depth + 1)
    return myTree

import pickle
# myTree = createTree([0,1], 0)
# print myTree

#inTree = open('myTree.dict')
#myTree = pickle.load(inTree)
#myTree['0, 1'][1] = '...'
#myTree['0, 1'][2] = '...'
#createPlot(myTree)
#inTree.close()

inTree = open('myTree.dict')
myTree = pickle.load(inTree)
createPlot(myTree)
inTree.close()

# outTree = open('myTree.dict', 'w')
# pickle.dump(myTree, outTree)
# outTree.close()
# createPlot(myTree)