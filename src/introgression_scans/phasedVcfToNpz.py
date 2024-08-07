import sys, random
import numpy as np

"""
args:
    vcfFileName -- the input vcf file
    species1ListFileName -- a file with a list of all individuals in the first population/species (names need to match identifiers in the vcf header)
    species2ListFileName -- a file with a list of all individuals in the second population/species (names need to match identifiers in the vcf header)
    maskedRefFileName -- a fasta file with the reference genome, and with all sites that we wish to ignored being masked as Ns
    arm -- the name of the chromosome/arm to be examined in this run (needs to match a corresponding entry in the fasta file specified by maskedRefFileName)
    npzFileName -- the output file to be fed to the predictor
"""
vcfFileName, species1ListFileName, species2ListFileName, maskedRefFileName, arm, npzFileName = sys.argv[1:]

def readSpeciesSampleList(sampleFileName):
    samps = {}
    with open(sampleFileName) as f:
        for line in f:
            samps[line.strip()] = 1
    return samps

def readFa(faFileName, upper=False):
    seqData = {}
    with open(faFileName) as faFile:
        reading = False
        for line in faFile:
            if line.startswith(">"):
                if reading:
                    if upper:
                        seqData[currChr] = seq.upper()
                    else:
                        seqData[currChr] = seq
                else:
                    reading = True
                currChr = line[1:].split()[0]
                seq = ""
            else:
                seq += line.strip()
    if upper:
        seqData[currChr] = seq.upper()
    else:
        seqData[currChr] = seq
    return seqData

def getDiploidHaps(genos, ref, alt, maskBase):
    allele1, allele2 = genos.split("|")
    #if allele1 == "." or maskBase == 'N':
    #    allele1 = -1
    #if allele2 == "." or maskBase == 'N':
    #    allele2 = -1
    #    allele1 = -1
    #if allele2 == "." or maskBase == 'N':
    #    allele2 = -1
    assert allele1 != "."
    assert allele2 != "."
    assert maskBase != 'N'
    return int(allele1), int(allele2)

def readSnpHapsFromPhasedSimSechVcf(vcfFileName, maskedRef, arm, spec1List, spec2List):
    simIndices, sechIndices = [], []
    simNames, sechNames = [], []
    snpGenosSim, snpGenosSech = {}, {}
    with open(vcfFileName) as vcfFile:
        for line in vcfFile:
            if line.startswith("#CHROM"):
                line = line.strip().split()
                for i in range(9, len(line)):
                    if line[i] in spec1List:
                        sechIndices.append(i)
                        sechNames.append(line[i]+".1")
                        sechNames.append(line[i]+".2")
                    elif line[i] in spec2List:
                        simIndices.append(i)
                        simNames.append(line[i]+".1")
                        simNames.append(line[i]+".2")
                header = line
            elif not line.startswith("#"):
                line = line.strip().split()
                c, pos, varId, ref, alt = line[:5]
                if c == arm:
                    assert ref in ['G','T','C','A'] and alt in ['G','T', 'C', 'A']
                    pos = int(pos)
                    snpGenosSim[(c, pos)] = {}
                    snpGenosSech[(c, pos)] = {}
                    for i in simIndices:
                        snpGenosSim[(c, pos)][header[i]+".1"], snpGenosSim[(c, pos)][header[i]+".2"] = getDiploidHaps(line[i], ref, alt, maskedData[c][pos-1])
                    for i in sechIndices:
                        snpGenosSech[(c, pos)][header[i]+".1"], snpGenosSech[(c, pos)][header[i]+".2"] = getDiploidHaps(line[i], ref, alt, maskedData[c][pos-1])
    return simNames, snpGenosSim, sechNames, snpGenosSech

def getMinMaj(genos1, genos2):
    allGenos = list(genos1.values()) + list(genos2.values())
    count0 = allGenos.count(0)
    count1 = allGenos.count(1)
    assert count0+count1 == len(allGenos)
    if count0 < count1:
        return [0, 1]
    elif count1 < count0:
        return [1, 0]
    elif count0 == count1:
        ls = [0,1]
        random.shuffle(ls)
        return ls

def recodeMinMaj(geno, minorAllele, majorAllele):
    if geno == minorAllele:
        return 1
    elif geno == majorAllele:
        return 0
    else:
        raise Exception

def writeNpzFile(headersSim, snpGenosSim, headersSech, snpGenosSech, npzFileName):
    simKeys = sorted(snpGenosSim.keys())
    sechKeys = sorted(snpGenosSech.keys())
    assert simKeys == sechKeys

    arms = {}
    simMatrix = []
    sechMatrix = []
    positions = []
    for c, pos in simKeys:
        minorAllele, majorAllele = getMinMaj(snpGenosSim[(c,pos)], snpGenosSech[(c,pos)])
        arms[c] = 1
        simVector = []
        for g in headersSim:
            simVector.append(recodeMinMaj(snpGenosSim[(c,pos)][g], minorAllele, majorAllele))

        sechVector = []
        for g in headersSech:
            sechVector.append(recodeMinMaj(snpGenosSech[(c,pos)][g], minorAllele, majorAllele))

        simMatrix.append(simVector)
        sechMatrix.append(sechVector)
        positions.append(pos)
    assert len(arms) == 1

    headersSim = np.array(headersSim)
    headersSech = np.array(headersSech)
    simMatrix = np.array(simMatrix, dtype='uint8')
    sechMatrix = np.array(sechMatrix, dtype='uint8')
    positions = np.array(positions, dtype='int')
    assert simMatrix.shape == (len(positions), len(headersSim))
    assert sechMatrix.shape == (len(positions), len(headersSech))

    np.savez_compressed(npzFileName, simHeader=headersSim, sechHeader=headersSech,
             simMatrix=simMatrix, sechMatrix=sechMatrix, positions=positions)

#read the list of sample names for each species
spec1Samps = readSpeciesSampleList(species1ListFileName)
spec2Samps = readSpeciesSampleList(species2ListFileName)

#read in the reference genome, which should have any sites that we want to ignore masked by 'N'
maskedData = readFa(maskedRefFileName, upper=True)

#read our vcf and extract the genotypes for each species
headersSim, snpGenosSim, headersSech, snpGenosSech = readSnpHapsFromPhasedSimSechVcf(vcfFileName, maskedData, arm, spec1Samps, spec2Samps)

#write the resulting header and genotype data as python/numpy objects in an npz file
writeNpzFile(headersSim, snpGenosSim, headersSech, snpGenosSech, npzFileName)
