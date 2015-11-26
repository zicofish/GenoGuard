'''
Created on 2014-2-9

@author: lenovo
'''
from __future__ import print_function
from collections import OrderedDict
from numpy import *

def alleleFrequencyProcess(alleleFreqsPath):
    alleleFreqs = {}
    afFile = open(alleleFreqsPath)
    header = afFile.readline()
    charSet = {"T", "A", "C", "G"}
    line = afFile.readline()
    while line != "":
        allele = line.split()
        line = afFile.readline()
        if float(allele[-6]) > float(allele[-3]):
            majorAlleleIndex = -7
            minorAlleleIndex = -4
        else:
            majorAlleleIndex = -4
            minorAlleleIndex = -7
        if allele[majorAlleleIndex] not in charSet or allele[minorAlleleIndex] not in charSet:
            continue
        alleleFreqs[int(allele[2])] = [allele[majorAlleleIndex], float(allele[majorAlleleIndex + 1]), allele[minorAlleleIndex], float(allele[minorAlleleIndex + 1])]
    afFile.close()
    return alleleFreqs

#range(100000000) needs 1.6G
#range(10000000) needs 165M
def ldProcess(ldPath, alleleFreqs):
    ld = []
    ldAlleles = set()
    ldFile = open(ldPath)
    firstAlleleIndex = 0
    secondAlleleIndex = 1
    DPrimeIndex = 5
    r2Index = 6
    line = ldFile.readline()
    while line != "":
        record = line.split()
        line = ldFile.readline()
        pos1 = int(record[firstAlleleIndex])
        pos2 = int(record[secondAlleleIndex])
        if (not alleleFreqs.has_key(pos1)) or (not alleleFreqs.has_key(pos2)):
            continue
        r2 = float(record[r2Index])
        p = [alleleFreqs[pos1][1], alleleFreqs[pos1][3]]
        q = [alleleFreqs[pos2][1], alleleFreqs[pos2][3]]
        D_abs = math.sqrt(r2*p[0]*p[1]*q[0]*q[1])
        
        DPrime = float(record[DPrimeIndex])
        
        positive = math.fabs(D_abs / min(p[0]*q[1], p[1]*q[0]) - DPrime)
        negative = math.fabs(D_abs / min(p[0]*q[0], p[1]*q[1]) - DPrime)
        if positive < negative:
            D = D_abs
        else:
            D = -D_abs
        ld.append([pos1, pos2, D])
        ldAlleles.add(pos1)
        ldAlleles.add(pos2)
    ldFile.close()
    return ld, ldAlleles

def genotypesProcess(genotypesPath, alleleFreqs, ldAlleles):
    genotypes = {}
    gtFile = open(genotypesPath)
    header = gtFile.readline().split()
    alleleValueIndex = 1
    posIndex = 3
    patientIndex = 0
    while not (header[patientIndex].startswith('NA')):
        patientIndex += 1
    line = gtFile.readline()
    while line != "":
        record = line.split()
        line = gtFile.readline()
        pos = int(record[posIndex])
        if not alleleFreqs.has_key(pos) or pos not in ldAlleles:
            continue
        alleleValues = record[alleleValueIndex].split('/')
        if alleleFreqs[pos][0] == alleleValues[1]:
            tmp = alleleValues[1]
            alleleValues[1] = alleleValues[0]
            alleleValues[0] = tmp
        for i in range(patientIndex, len(header)):
            if record[i] == 'NN':
                tmp0 = alleleFreqs[pos][0] if random.random() < alleleFreqs[pos][1] else alleleFreqs[pos][2]
                tmp1 = alleleFreqs[pos][0] if random.random() < alleleFreqs[pos][1] else alleleFreqs[pos][2]
                record[i] = tmp0+tmp1
        convertToNum = {alleleValues[0]*2:0, alleleValues[0]+alleleValues[1]:1, alleleValues[1]+alleleValues[0]:1, alleleValues[1]*2:2}
        try:
            genotypes[pos] = [convertToNum[record[i]] for i in range(patientIndex, len(header))]
        except KeyError:
            print(pos)
            return
    gtFile.close()
    return genotypes
    
def saveNewData(genotypes, alleleFreqs, ld, newGenotypesPath, newAlleleFreqsPath, newLdPath):
    newGtFile = open(newGenotypesPath, 'w')
    newAfFile = open(newAlleleFreqsPath, 'w')
    newLdFile = open(newLdPath, 'w')
    
    newGtFile.write("\n".join(["\t".join([str(i) for i in [item[0]] + item[1]]) for item in sorted(genotypes.items(), key = lambda item: item[0])]))
    genotypes = {}
    newAfFile.write("\n".join(["\t".join([str(i) for i in [item[0]] + item[1]]) for item in sorted(alleleFreqs.items(), key = lambda item: item[0])]))
    alleleFreqs = {}
    for ld_item in ld:
        newLdFile.write("\t".join([str(i) for i in ld_item]) + "\n")
    ld = []
    newGtFile.close()
    newAfFile.close()
    newLdFile.close()

def processAllChroms():
    for i in range(10, 11):
        alleleFreqs = {}
        ld = []
        ldAlleles = set()
        genotypes = {}
        alleleFreqs = alleleFrequencyProcess("./hapmap/chr%s/allele_freqs_chr%s_CEU_phase3.2_nr.b36_fwd.txt" % (i, i))
        ld, ldAlleles = ldProcess("./hapmap/chr%s/ld_chr%s_CEU.txt" % (i, i), alleleFreqs)
        genotypes = genotypesProcess("./hapmap/chr%s/genotypes_chr%s_CEU_phase3.2_consensus.b36_fwd.txt" % (i, i), alleleFreqs, ldAlleles)
        
        saveNewData(genotypes, alleleFreqs, ld, "./hapmap/chr%s/big_genotypes_chr%s_CEU.txt" % (i, i), "./hapmap/chr%s/tmp_allele_freqs_chr%s_CEU.txt" % (i, i), "./hapmap/chr%s/tmp_ld_chr%s_CEU.txt" % (i, i))


if __name__ == '__main__':
    processAllChroms()
#     alleleFreqs = alleleFrequencyProcess("./hapmap/allele_freqs_chr22_CEU_phase3.2_nr.b36_fwd.txt")
#     ld, ldAlleles = ldProcess("./hapmap/ld_chr22_CEU.txt", alleleFreqs)
#     genotypes = genotypesProcess("./hapmap/genotypes_chr22_CEU_phase3.2_consensus.b36_fwd.txt", alleleFreqs, ldAlleles)
#     saveNewData(genotypes, alleleFreqs, ld, "./hapmap/new_genotypes_chr22_CEU.txt", "./hapmap/new_allele_freqs_chr22_CEU.txt", "./hapmap/new_ld_chr22_CEU.txt")
