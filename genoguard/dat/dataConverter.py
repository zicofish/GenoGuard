'''
Created on Jun 16, 2014

@author: zhihuang
'''
top100SNPs = ['rs194014', 'rs2297679', 'rs7013027', 'rs733290', 'rs1844396',
'rs1356733', 'rs975664', 'rs2822987', 'rs10134505', 'rs7539636',
'rs4751629', 'rs152280', 'rs2324909', 'rs6050374', 'rs4785919',
'rs1281340', 'rs156984', 'rs3943733', 'rs1512744', 'rs641902',
'rs11164559', 'rs4242165', 'rs1463652', 'rs12670839', 'rs2589655',
'rs3128917', 'rs2303942', 'rs10505909', 'rs441517', 'rs7358335',
'rs10188217', 'rs2427622', 'rs10746709', 'rs2276070', 'rs6681880',
'rs10872465', 'rs2424303', 'rs3873386', 'rs4870268', 'rs867926',
'rs1467044', 'rs7656234', 'rs3105439', 'rs2417680', 'rs4577244',
'rs2172159', 'rs746124', 'rs9314379', 'rs10519410', 'rs11096674',
'rs6583859', 'rs6741951', 'rs10497705', 'rs39108', 'rs4973190',
'rs9960403', 'rs10510633', 'rs3006198', 'rs6958889', 'rs4269797',
'rs2441727', 'rs10189760', 'rs16913731', 'rs7732983', 'rs181247',
'rs12777190', 'rs4800827', 'rs10850434', 'rs10515467', 'rs11195943',
'rs4706511', 'rs1391451', 'rs11811536', 'rs12576387', 'rs7660783',
'rs12759306', 'rs10493856', 'rs17494100', 'rs7404238', 'rs1322755',
'rs11097457', 'rs2256965', 'rs1801249', 'rs149332', 'rs7580688',
'rs12594390', 'rs4280690', 'rs4648527', 'rs12347149', 'rs1906990',
'rs7646054', 'rs2102272', 'rs9812253', 'rs4741546', 'rs2918285',
'rs6012906', 'rs3816969', 'rs6496173', 'rs4921809', 'rs250850']

def mapFileTrunc(oldMapFileName, newMapFileName):
    mapFile = open(oldMapFileName)
    newMapFile = open(newMapFileName, 'w')
    for line in mapFile.readlines():
        attrArray = line.split()
        if len(attrArray) >= 4 and attrArray[1] in top100SNPs:
            newMapFile.write("\t".join(attrArray) + "\n")
    mapFile.close()
    newMapFile.close()
    
def pedFileTrunc(oldPedFileName, newPedFileName, oldMapFileName):
    pedFile = open(oldPedFileName, buffering = 1)
    newPedFile = open(newPedFileName, mode = 'w', buffering = 1)
    mapFile = open(oldMapFileName)
    chosenSNPLines = []
    lineNum = -1
    for line in mapFile.readlines():
        lineNum += 1
        attrArray = line.split()
        if len(attrArray) >= 4 and attrArray[1] in top100SNPs:
            chosenSNPLines.append(lineNum)
    line = pedFile.readline()
    while line != '':
        attrArray = line.split()
        newAttrArray = attrArray[:6]
        for lineNum in chosenSNPLines:
            newAttrArray.append(attrArray[lineNum*2 + 6])
            newAttrArray.append(attrArray[lineNum*2 + 7])
        newPedFile.write('\t'.join(newAttrArray) + "\n")
        line = pedFile.readline()
    mapFile.close()
    newPedFile.close()
    pedFile.close()
   
def getMajorAlleles(SNPRefListFileName):
    SNPRefList = open(SNPRefListFileName)
    majorAlleles = [line.split()[2] for line in SNPRefList.readlines()]
    return majorAlleles

def getAncestryMap(ancestryFileName):
    ancestryFile = open(ancestryFileName)
    ancestryMap = {}
    for line in ancestryFile.readlines():
        attrArray = line.split()
        ancestryMap[attrArray[1]] = attrArray[6]
    return ancestryMap
    
def convertToPCATrainingSet(pedFileName, mapFileName, ancestryFileName, 
                            SNPRefListFileName, outFileName):
    majorAlleles = getMajorAlleles(SNPRefListFileName)
    ancestryMap = getAncestryMap(ancestryFileName)
    pedFile = open(pedFileName)
    mapFile = open(mapFileName)
    trainingSet = []
    for line in pedFile.readlines():
        attrArray = line.split()
        dataPoint = [' '] * (len(majorAlleles) + 2)
        dataPoint[0] = attrArray[1]
        dataPoint[1] = ancestryMap[attrArray[1]]
        dataPoint[2:] = [str(2 - (attrArray[2*i+6]+attrArray[2*i+7]).count(majorAlleles[i])) 
                         for i in range(len(majorAlleles))]
        trainingSet.append(dataPoint)
    outFile = open(outFileName, 'w')
    outFile.write('\n'.join(['\t'.join(point) for point in trainingSet]))
    outFile.close()
    pedFile.close()
    mapFile.close()

def SNPList(mapFileName, outFileName):
    mapFile = open(mapFileName)
    outFile = open(outFileName, 'w')
    snps = [line.split()[1] for line in mapFile.readlines()]
    outFile.write('\n'.join(snps))
    outFile.close()
    mapFile.close()
    
def alleleFreqCompile(freqFolderName, mapFileName, ancestry):
    import os
    files = [f for f in os.listdir(freqFolderName) 
             if f.find(ancestry) >= 0]
    mapFile = open(mapFileName)
    plinkMap = mapFile.readlines()
    snps = [line.split()[1] for line in plinkMap]
    outData = []
    for f in files:
        freqFile = open(freqFolderName + f)
        freqFile.readline()  #header
        for line in freqFile.readlines():
            attrArray = line.split()
            if attrArray[0] not in snps:
                continue
#             if float(attrArray[-6]) > float(attrArray[-3]):
#                 majorAlleleIndex = -7
#                 minorAlleleIndex = -4
#             else:
#                 majorAlleleIndex = -4
#                 minorAlleleIndex = -7
            outData.append([attrArray[1], attrArray[0], attrArray[2],
                            attrArray[-7],
                            attrArray[-6],
                            attrArray[-4],
                            attrArray[-3]])
        freqFile.close()
    outData = sorted(outData, key=lambda x: snps.index(x[1]))
    outFile = open(ancestry + '_freq.txt', 'w')
    outFile.write('\n'.join(['\t'.join(line) for line in outData]))
    outFile.close()
    mapFile.close()

def genSNPRefList(freqFileNameList, outFileName):
    outList = []
    SNPRefLists = []
    for freqFileName in freqFileNameList:
        freqFile = open(freqFileName)
        SNPRefLists.append([line.split() for line in freqFile.readlines()])
        freqFile.close()
    for i in range(len(SNPRefLists[0])):
        majors = [refList[i][3] for refList in SNPRefLists]
        if majors.count(majors[0]) == len(majors):
            outList.append([SNPRefLists[0][i][0], SNPRefLists[0][i][1], 
                            SNPRefLists[0][i][3], SNPRefLists[0][i][5]])
        else:
            raise ValueError("Non consensus major at: " + SNPRefLists[0][1])
    outFile = open(outFileName, 'w')
    outFile.write('\n'.join(['\t'.join(SNP) for SNP in outList]))
    outFile.close()
        
def genChrom22Map(oldMapFileName, newMapFileName):
    mapFile = open(oldMapFileName)
    newMapFile = open(newMapFileName, 'w')
    for line in mapFile.readlines():
        attrArray = line.split()
        if len(attrArray) >= 4 and attrArray[0] == '22':
            newMapFile.write("\t".join(attrArray) + "\n")
    mapFile.close()
    newMapFile.close()
            
def genChrom22Ped(oldPedFileName, newPedFileName, oldMapFileName):
    pedFile = open(oldPedFileName, buffering = 1)
    newPedFile = open(newPedFileName, mode = 'w', buffering = 1)
    mapFile = open(oldMapFileName)
    chosenSNPLines = []
    lineNum = -1
    for line in mapFile.readlines():
        lineNum += 1
        attrArray = line.split()
        if len(attrArray) >= 4 and attrArray[0] == '22':
            chosenSNPLines.append(lineNum)
    line = pedFile.readline()
    while line != '':
        attrArray = line.split()
        newAttrArray = attrArray[:6]
        for lineNum in chosenSNPLines:
            newAttrArray.append(attrArray[lineNum*2 + 6])
            newAttrArray.append(attrArray[lineNum*2 + 7])
        newPedFile.write('\t'.join(newAttrArray) + "\n")
        line = pedFile.readline()
    mapFile.close()
    newPedFile.close()
    pedFile.close()

'''
Filter the snps in the allele frequency file so that it has the same SNPList as the genotype file
'''
def pruneFreqList(oldFreqFileName, newFreqFileName, datFileName):
    datFile = open(datFileName)
    freqFile = open(oldFreqFileName)
    outFile = open(newFreqFileName, 'w')
    SNPList = [line.split()[0] for line in datFile.readlines()]
    freqs = []
    j = 0
    for line in freqFile.readlines():
        attrArray = line.split()
        while j < len(SNPList) and int(attrArray[0]) > int(SNPList[j]):
            j += 1
        if j >= len(SNPList):
            break
        elif int(attrArray[0]) == int(SNPList[j]):
            freqs.append(attrArray)
    outFile.write('\n'.join(['\t'.join(line) for line in freqs]))
    outFile.close()
    freqFile.close()
    datFile.close()
    
'''
Filter LDs so that it has the same SNPList as the genotype file, and for each snp, retain at most one LD with previous snps
'''
def pruneLD(oldLDFileName, newLDFileName, datFileName):
    datFile = open(datFileName)
    LDFile = open(oldLDFileName)
    outFile = open(newLDFileName, 'w')
    SNPList = [line.split()[0] for line in datFile.readlines()]
    prefixLD = {}
    for line in LDFile.readlines():
        attrArray = line.split()
        # Only retains the previous closest LD for each snp
        prefixLD[attrArray[1]] = [attrArray[0], attrArray[2]]
    prefixLDArr = []
    for idx in range(len(SNPList)):
        if prefixLD.has_key(SNPList[idx]):
            try:
                idx_tmp = SNPList.index(prefixLD[SNPList[idx]][0])
                prefixLDArr.append([SNPList[idx], str(idx_tmp), prefixLD[SNPList[idx]][1]])
            except ValueError:
                # The previous snp does not exist in the genotype SNPList
                prefixLDArr.append([SNPList[idx]])
        else:
            # This snp has no LD with any previous snp
            prefixLDArr.append([SNPList[idx]])
    outFile.write('\n'.join(['\t'.join(line) for line in prefixLDArr]))
    outFile.close()
    LDFile.close()
    datFile.close()

'''

'''    
def processRecombData(AFFileName, recombFileName, legendFileName, haplotypeFileName, newRecombFileName, newHaplotypeFileName):
    AFFile = open(AFFileName)
    AF = []
    SNPRefList = []
    for line in AFFile.readlines():
        attrArray = line.split()
        AF.append([attrArray[1], attrArray[3]])  # I only need the letters for major allele and minor allele, but not their frequencies
        SNPRefList.append(attrArray[0])
    AFFile.close()

    recombFile = open(recombFileName)
    recombFile.readline()  # header
    recomb = map(lambda u: u.split(), recombFile.readlines())
    recombFile.close()

    legendFile = open(legendFileName)
    legendFile.readline()  # header
    legend = map(lambda u: u.split()[1:], legendFile.readlines())  # I don't need the rsID in the first column
    legendFile.close()
    
    haplotypeFile = open(haplotypeFileName)
    haplotype = map(lambda u: map(int, u.split()), haplotypeFile.readlines())
    haplotypeFile.close()
    
    # Validation phase... Make sure that every SNP in the SNPRefList exists also in legend
    print "Start validation........................."
    legendIdx = 0
    for snp in SNPRefList:
        flag = False
        for i in range(legendIdx, len(legend)):
            if legend[i][0] == snp:
                flag = True
                legendIdx = i + 1
                break
        if not flag:
            raise ValueError("snp %s not found in legend" % (snp))
    print "Validation succeed..."
    
    # Filter and transform haplotype dat. For each snp, if the two alleles in the legend file has different order from those 
    # in the allele frequency file, we should switch them in the legend file and changes 0 to 1 and 1 to 0 in the haplotype file
    print "Start filtering and transforming haplotype dat.................."
    newHaplotype = []
    legendIdx = 0
    for i in range(len(SNPRefList)):
        for j in xrange(legendIdx, len(legend)):
            if legend[j][0] == SNPRefList[i]:
                if legend[j][1] == AF[i][0]:  # The two alleles have the same order in the legend file as in the allele frequency file
                    newHaplotype.append([SNPRefList[i]] + haplotype[j])
                else:  # Different order
                    print "snp %s has a different allele order" % (SNPRefList[i])
                    newHaplotype.append([SNPRefList[i]] + map(lambda u: 1-u, haplotype[j]))  # 0 to 1 and 1 to 0
                legendIdx = j+1
                break
    print "Haplotype transformation finished........."
    
    # Preprocess the recombination rates file, so that it is convenient for our use.
    # Namely, we would construct a new recombination rates file, where each row is of the form "position genetic_distance"
    # The list of position is the same as in allele frequency file.
    # The genetic distance is between this snp and the previous snp (recombination_rate * physical distance)
    print "Start processing recombination rates..................."
    newRecomb = []
    newRecomb.append([SNPRefList[0], 0])
    recombIdx = -1
    while int(recomb[recombIdx+1][0]) < int(SNPRefList[0]):
        recombIdx += 1
    i = 1
    while i < len(SNPRefList):
        pos = int(SNPRefList[i-1])
        distance = 0
        while (recombIdx+1) < len(recomb) and int(recomb[recombIdx + 1][0]) < int(SNPRefList[i]):
            if float(recomb[recombIdx][1]) < 0:
                recomb[recombIdx][1] = 0
            distance += float(recomb[recombIdx][1]) * (int(recomb[recombIdx + 1][0]) - pos) / 1000000
            pos = int(recomb[recombIdx + 1][0])
            recombIdx += 1
        if float(recomb[recombIdx][1]) < 0:
            recomb[recombIdx][1] = 0
        distance += float(recomb[recombIdx][1]) * (int(SNPRefList[i]) - pos) / 1000000
        newRecomb.append([SNPRefList[i], distance])
        i += 1
    print "Recombination rates processing finished..............."
    
    print "Writing out dat.............."
    newRecombFile = open(newRecombFileName, 'w')
    newRecombFile.write('\n'.join(map(lambda u: '\t'.join(map(str, u)), newRecomb)))
    newRecombFile.close()
    newHaplotypeFile = open(newHaplotypeFileName, 'w')
    newHaplotypeFile.write('\n'.join(map(lambda u: '\t'.join(map(str, u)), newHaplotype)))
    newHaplotypeFile.close()
    print "Done"
    
def trunkToSmall(genotypeFileName, newGenotypeFileName, AFFileName, newAFFileName, LDFileName, newLDFileName, recombFileName, newRecombFileName, haplotypeFileName, newHaplotypeFileName):
    genotypeFile = open(genotypeFileName)
    AFFile = open(AFFileName)
    LDFile = open(LDFileName)
    recombFile = open(recombFileName)
    haplotypeFile = open(haplotypeFileName)
    smallGenotypeFile = open(newGenotypeFileName, 'w')
    smallAFFile = open(newAFFileName, 'w')
    smallLDFile = open(newLDFileName, 'w')
    smallRecombFile = open(newRecombFileName, 'w')
    smallHaplotypeFile = open(newHaplotypeFileName, 'w')
    smallGenotypeFile.write(''.join(genotypeFile.readlines()[:1000]))
    smallAFFile.write(''.join(AFFile.readlines()[:1000]))
    smallLDFile.write(''.join(LDFile.readlines()[:1000]))
    smallRecombFile.write(''.join(recombFile.readlines()[:1000]))
    smallHaplotypeFile.write(''.join(haplotypeFile.readlines()[:1000]))
    smallGenotypeFile.close()
    smallAFFile.close()
    smallLDFile.close()
    smallRecombFile.close()
    smallHaplotypeFile.close()
    
def pruneFileToolkit():
    import os
    for i in range(5, 23):
        os.remove('./hapmap/chr%s/allele_freqs_chr%s_CEU_phase3.2_nr.b36_fwd.txt' % (i, i))
        os.remove('./hapmap/chr%s/CEU.chr%s.hap' % (i, i))
        os.remove('./hapmap/chr%s/genetic_map_chr%s_combined_b36.txt' % (i, i))
        os.remove('./hapmap/chr%s/genotypes_chr%s_CEU_phase3.2_consensus.b36_fwd.txt' % (i, i))
        os.remove('./hapmap/chr%s/hapmap3.r2.b36.chr%s.legend' % (i, i))
        os.remove('./hapmap/chr%s/ld_chr%s_CEU.txt' % (i, i))
        os.remove('./hapmap/chr%s/new_allele_freqs_chr%s_CEU.txt' % (i, i))
        os.remove('./hapmap/chr%s/new_ld_chr%s_CEU.txt' % (i, i))
        
        os.rename('./hapmap/chr%s/hg_allele_freqs_chr%s_CEU.txt' % (i, i), './hapmap/chr%s/big_allele_freqs_chr%s_CEU.txt' % (i, i))
        os.rename('./hapmap/chr%s/hg_ld_chr%s_CEU.txt' % (i, i), './hapmap/chr%s/big_ld_chr%s_CEU.txt' % (i, i))
        os.rename('./hapmap/chr%s/new_CEU.chr%s.hap' % (i, i), './hapmap/chr%s/big_CEU.chr%s.hap' % (i, i))
        os.rename('./hapmap/chr%s/new_genetic_map_chr%s_combined_b36.txt' % (i, i), './hapmap/chr%s/big_genetic_map_chr%s_combined_b36.txt' % (i, i))
        os.rename('./hapmap/chr%s/new_genotypes_chr%s_CEU.txt' % (i, i), './hapmap/chr%s/big_genotypes_chr%s_CEU.txt' % (i, i))
        
def computeAllChromsIntSeqs():
    import gmpy2
    import numpy as np
    import pickle
    for i in range(1, 23):
        dataset = []
        datFile = open('./hapmap/chr%s/big_genotypes_chr%s_CEU.txt' % (i, i))
        for line in datFile.readlines():
            attrArray = line.split()
            dataset.append(attrArray[1:])
        dataset = np.array(dataset, dtype=int)
        intSeqs = []
        for j in range(np.shape(dataset)[1]):
            intSeqs.append(gmpy2.mpz(''.join(map(str, list(dataset[:, j]))), 3))
        intSeqs = sorted(intSeqs)  # the intseqs are sorted, so they don't correspond to the original patients in order
        
        intSeqsFile = open('./hapmap/chr%s/big_int_seqs_chr%s_CEU.pickle' % (i, i), 'w')
        pickle.dump(intSeqs, intSeqsFile)
        intSeqsFile.close()
        

def convertAllChroms():
    for i in range(10, 11):
        pruneFreqList("./hapmap/chr%s/tmp_allele_freqs_chr%s_CEU.txt" % (i, i),
                      "./hapmap/chr%s/big_allele_freqs_chr%s_CEU.txt" % (i, i),
                      "./hapmap/chr%s/big_genotypes_chr%s_CEU.txt" % (i, i))
        pruneLD("./hapmap/chr%s/tmp_ld_chr%s_CEU.txt" % (i, i),
                "./hapmap/chr%s/big_ld_chr%s_CEU.txt" % (i, i),
                "./hapmap/chr%s/big_genotypes_chr%s_CEU.txt" % (i, i))
        processRecombData('./hapmap/chr%s/big_allele_freqs_chr%s_CEU.txt' % (i, i),
                          './hapmap/chr%s/genetic_map_chr%s_combined_b36.txt' % (i, i),
                          './hapmap/chr%s/hapmap3.r2.b36.chr%s.legend' % (i, i),
                          './hapmap/chr%s/CEU.chr%s.hap' % (i, i),
                          './hapmap/chr%s/big_genetic_map_chr%s_combined_b36.txt' % (i, i),
                          './hapmap/chr%s/big_CEU.chr%s.hap' % (i, i))
        trunkToSmall('./hapmap/chr%s/big_genotypes_chr%s_CEU.txt' % (i, i),
                     './hapmap/chr%s/small_genotypes_chr%s_CEU.txt' % (i, i),
                     './hapmap/chr%s/bigg_allele_freqs_chr%s_CEU.txt' % (i, i),
                     './hapmap/chr%s/small_allele_freqs_chr%s_CEU.txt' % (i, i),
                     './hapmap/chr%s/bigg_ld_chr%s_CEU.txt' % (i, i),
                     './hapmap/chr%s/small_ld_chr%s_CEU.txt' % (i, i),
                     './hapmap/chr%s/big_genetic_map_chr%s_combined_b36.txt' % (i, i),
                     './hapmap/chr%s/small_genetic_map_chr%s_combined_b36.txt' % (i, i),
                     './hapmap/chr%s/big_CEU.chr%s.hap' % (i, i),
                     './hapmap/chr%s/small_CEU.chr%s.hap' % (i, i))
    
def transposeData(inFileName, outFileName):
    import numpy
    inFile = open(inFileName)
    datArr = map(lambda u: map(int, u.split()), inFile.readlines())
    inFile.close()
    datMat = numpy.mat(datArr).transpose().getA()
    outFile = open(outFileName, 'w')
    outFile.write("\n".join(map(lambda u: " ".join(map(str, u)), datMat)))
    outFile.close()
    
'''
Randomly hide the specified percentage of SNVs in the input genotype file
'''
def hideSNVs(kwargs):
    import random
    genotypeFile = open(kwargs['genotypeFileName'])
    genotype = map(lambda u: u.split()[kwargs['startIdx']:kwargs['endIdx']], genotypeFile.readlines())
    hiddenGenotype = []
    for i in range(len(genotype)): # length of genotype
        tmp = []
        for j in range(len(genotype[0])): # number of patients
            r = random.random()
            if r < float(kwargs['percent']):
                tmp.append('H')
            else:
                tmp.append(genotype[i][j])
        hiddenGenotype.append(tmp)
    outFile = open(kwargs['outFileName'], 'w')
    outFile.write('\n'.join(map(lambda u: '\t'.join(u), hiddenGenotype)))
    outFile.close()
    
def alignHaploToGeno(kwargs):
    orgGenotypeFile = open(kwargs['orgGneotypeFileName'])
    header = orgGenotypeFile.readline().split()
    orgGenotypeFile.close()
    patientIndex = 0
    while not (header[patientIndex].startswith('NA')):
        patientIndex += 1
    pOrders = {}
    for i in range(patientIndex, len(header)):
        pOrders[header[i]] = i
        
    def m_cmp(x, y):
        xo = pOrders[x]
        yo = pOrders[y]
        if xo < yo:
            return -1
        if xo == yo:
            return 0
        return 1
    
    afFile = open(kwargs['afFileName'])
    af = {}
    snpList = []
    haploMat = []
    for line in afFile.readlines():
        strs = line.split()
        snpList.append(strs[0]) 
        af[strs[0]] = [strs[1], strs[3]]
        haploMat.append([])
    afFile.close()
    
    haploList = []
    for haploFileName in kwargs['haploFileNames']:
        haploFile = open(haploFileName)
        header = haploFile.readline().split()
        haploList.extend(header[2:])
        
        i = 0
        for line in haploFile.readlines():
            strs = line.split()
            if not af.has_key(strs[1]):
                continue
            haploMat[i].extend(map(lambda u: 0 if u == af[strs[1]][0] else 1, strs[2:]))
            i += 1
    import numpy
    haploMat = numpy.array(haploMat)
    haploDict = {}
    for i in range(len(haploList)):
        haploDict[haploList[i]] = list(haploMat[:, i])
    alignedHaplo = numpy.array(map(lambda u: u[1], sorted(haploDict.items(), cmp = lambda x, y: m_cmp(x[0].split('_')[0], y[0].split('_')[0])))).transpose()
    
    outFile = open(kwargs['outFileName'], 'w')
    outFile.write('Pos\t' + '\t'.join(sorted(haploList, cmp = lambda x, y: m_cmp(x.split('_')[0], y.split('_')[0]))) + '\n')
    outFile.write('\n'.join(map(lambda u: u[0] + '\t' + '\t'.join(map(str, u[1])), zip(snpList, alignedHaplo))))
    outFile.close()
    
if __name__ == '__main__':
    for i in range(1, 21):
        print "Hiding ", i*0.05*100, "% of SNVs and producing the result file \"hiddenSNVs_s0_e165_chr22 ", i*0.05, ".txt\""
        hideSNVs({'genotypeFileName':'../hapmap/chr22/small_genotypes_chr22_CEU.txt',
                  'startIdx':0,
                  'endIdx':165,
                  'percent':i*0.05,
                  'outFileName':'../hapmap/chr22/hiddenSNVs/hiddenSNVs_s0_e165_chr22_'+str(i*0.05)+'.txt'})
#     alignHaploToGeno({'orgGneotypeFileName':'./hapmap/chr22/genotypes_chr22_CEU_phase3.2_consensus.b36_fwd.txt',
#                       'afFileName':'./hapmap/chr22/big_allele_freqs_chr22_CEU.txt',
#                       'haploFileNames':['./hapmap/chr22/hapmap3_r2_b36_fwd.consensus.qc.poly.chr22_ceu.D.phased', 
#                                         './hapmap/chr22/hapmap3_r2_b36_fwd.consensus.qc.poly.chr22_ceu.phased',
#                                         './hapmap/chr22/hapmap3_r2_b36_fwd.consensus.qc.poly.chr22_ceu.unr.phased'],
#                       'outFileName':'./hapmap/chr22/big_aligned_CEU.chr22.hap'})
#     for i in range(1, 21):
#         transposeData("./merge-output/hapmap/chr22/infer_test/predict_recomb_chr22_CEU_%s.txt" % (i*0.05), 
#                       "./merge-output/hapmap/chr22/infer_test/PredictedSNVs_chr22_CEU_%s.txt" % (i*0.05))

#     transposeData("./hapmap/chr22/big_recomb_random_chr22_CEU.txt", 
#                   "./hapmap/chr22/synthetic_recomb_chr22_CEU.txt")

#     for i in range(10):
#         transposeData("./merge-output/hapmap/chr22/infer_test/predict_recomb_chr22_CEU_10percent_%s.txt" % (i), 
#                       "./merge-output/hapmap/chr22/infer_test/PredictedSNVs10percent_chr22_%s.txt" % (i))
#     for i in range(10):
#         transposeData("./merge-output/hapmap/chr22/infer_test/predict_recomb_chr22_CEU_40percent_%s.txt" % (i), 
#                       "./merge-output/hapmap/chr22/infer_test/PredictedSNVs40percent_chr22_%s.txt" % (i))

#     for i in range(10):
#         transposeData("./merge-output/hapmap/chr22/infer_training/predict_recomb_chr22_CEU_10percent_%s.txt" % (i), 
#                       "./merge-output/hapmap/chr22/infer_training/PredictedSNVs10percent_chr22_%s.txt" % (i))
#     for i in range(10):
#         transposeData("./merge-output/hapmap/chr22/infer_training/predict_recomb_chr22_CEU_40percent_%s.txt" % (i), 
#                       "./merge-output/hapmap/chr22/infer_training/PredictedSNVs40percent_chr22_%s.txt" % (i))
    
#     convertAllChroms()
#     pruneFileToolkit()
#     computeAllChromsIntSeqs()
    
# processRecombData()
    
# trunkToSmall('./hapmap/new_genotypes_chr22_CEU.txt',
#              './hapmap/hg_allele_freqs_chr22_CEU.txt',
#              './hapmap/hg_ld_chr22_CEU.txt')
#     pruneFreqList("./hapmap/new_allele_freqs_chr22_CEU.txt",
#                     "./hapmap/hg_allele_freqs_chr22_CEU.txt",
#                     "./hapmap/new_genotypes_chr22_CEU.txt")
# pruneLD("./hapmap/new_ld_chr22_CEU.txt",
#         "./hapmap/hg_ld_chr22_CEU.txt",
#         "./hapmap/new_genotypes_chr22_CEU.txt")

# genChrom22Map("C:/Users/zhihuang/Desktop/dat/hapmap/phase3/hapmap3_r2_b36_fwd.consensus.qc.poly.map",
#               "chrom22.map")

# genChrom22Ped("C:/Users/zhihuang/Desktop/dat/hapmap/phase3/hapmap3_r2_b36_fwd.consensus.qc.poly.ped",
#                  "chrom22.ped",
#                  "C:/Users/zhihuang/Desktop/dat/hapmap/phase3/hapmap3_r2_b36_fwd.consensus.qc.poly.map")
    
# mapFileConverter("C:/Users/zhihuang/Desktop/dat/hapmap/phase3/hapmap3_r2_b36_fwd.consensus.qc.poly.map", 
#                  "C:/Users/zhihuang/Desktop/dat/hapmap/phase3/top100SNPs.map")
# pedFileConverter("C:/Users/zhihuang/Desktop/dat/hapmap/phase3/hapmap3_r2_b36_fwd.consensus.qc.poly.ped",
#                  "C:/Users/zhihuang/Desktop/dat/hapmap/phase3/top100SNPs.ped",
#                  "C:/Users/zhihuang/Desktop/dat/hapmap/phase3/hapmap3_r2_b36_fwd.consensus.qc.poly.map")
# print getMajorAlleles("SNPRefList.txt")
# convertToPCATrainingSet("top100SNPs.ped",
#                        "top100SNPs.map",
#                        "relationships_w_pops_121708.txt",
#                        "SNPRefList.txt",
#                        "trainingSet.txt")
# SNPList("C:/Users/zhihuang/Desktop/dat/hapmap/phase3/top100SNPs.map",
#         "C:/Users/zhihuang/Desktop/dat/hapmap/phase3/query_snps.txt")
# alleleFreqCompile("C:/Users/zhihuang/Desktop/dat/hapmap/phase3/freq/",
#                   "top100SNPs.map",
#                   'CEU')
# genSNPRefList(['CEU_freq.txt', 'ASW_freq.txt', 'CHB_freq.txt'],
#               'SNPRefList.txt')
