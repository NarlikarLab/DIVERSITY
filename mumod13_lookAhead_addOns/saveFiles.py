import weblogoMod.weblogolib as wl
import numpy as np
import pickle
import os
import re
import plotFigures
#np.set_printoptions(threshold=np.inf)

def readData(dataFile):
    data = []
    tmpData = ""
    with open(dataFile) as infile:
        for line in infile:
            if line[0] == '>':
                if tmpData != "": data.append(tmpData)
                tmpData = ""
            else:
                tmpData = tmpData + line.strip()
    if tmpData != "": data.append(tmpData)
    return data

def getN(c):
    if c == 'N': return c
    elif c == 'A': return 'T'
    elif c == 'a': return 't'
    elif c == 'C': return 'G'
    elif c == 'c': return 'g'
    elif c == 'G': return 'C'
    elif c == 'g': return 'c'
    elif c == 'T': return 'A'
    else: return 'a'

def nToNum(c):
    if c == 'A' or c == 'a': return 0
    elif c == 'C' or c == 'c': return 1
    elif c == 'G' or c == 'g': return 2
    else: return 3

def reverseSequence(sequence):
    return ''.join(map(lambda x: getN(x), list(sequence))[::-1])

def createLogo(sequences, filename):
    if sequences == []: return
    seqs = wl.read_seq_data(sequences)
    data = wl.LogoData.from_seqs(seqs)
    options = wl.LogoOptions()
    options.logo_title = str(len(sequences)) + " sequences"
    options.size = "large"
    options.color_scheme = wl.colorscheme.monochrome
    formt = wl.LogoFormat(data, options)
    fout = open(filename, "w")
    wl.png_formatter(data, formt, fout)
    fout.close()

def getPSSM(motif):
    n = len(motif)
    if n == 0: return
    mWidth = len(motif[0])
    motifArr = np.array(map(lambda x: map(lambda y: nToNum(y), list(x)), motif)).transpose()
    motifArr = map(lambda x: np.bincount(x, minlength = 4), motifArr)
    motifArr = map(lambda x: map(lambda y: y/(1.0 * n), x), motifArr)
    motifArr = np.array(motifArr).transpose()
    return motifArr

def savePSSMMode(motifs, n, mode, filename, likelihood):
    motifArrays = map(lambda x: getPSSM(x), motifs)
    arr = ["A", "C", "G", "T"]
    f = open(filename, "w")
    f.write("\nNumber of sequences: " + str(n) + "\nNumber of modes: " + str(mode) + "\nLikelihood: " + str(likelihood) + "\n\n")
    for i in range(mode):
        f.write("\nMOTIF: " + str(i) + " (" + str(len(motifs[i])) + " sequences, width ")
        if motifArrays[i] is None: 
            f.write("0)\n\n")
            continue
        mLen = len(motifArrays[i][0])
        f.write(str(mLen) + ")\n\n")
        for j in range(4):
            f.write(arr[j] + " [")
            for k in range(mLen-1):
                f.write(str(motifArrays[i][j][k]) + "\t")
            f.write(str(motifArrays[i][j][mLen - 1]) + "]\n")
        f.write("\n\n")
    #     f.write("A [")
    #     for j in range(len(motifArrays[i])):
    #         f.write(str(j).zfill(2) + "\t" + str(motifArrays[i][j][0]) + "\t" + str(motifArrays[i][j][1]) + "\t" + str(motifArrays[i][j][2]) + "\t" + str(motifArrays[i][j][3]) + "\n")


    #     f.write("\n\n")

# def getPSSM(motif):
#     n = len(motif)
#     if n == 0: return
#     mWidth = len(motif[0])
#     motifArr = np.array(map(lambda x: map(lambda y: nToNum(y), list(x)), motif)).transpose()
#     motifArr = map(lambda x: np.bincount(x, minlength = 4), motifArr)
#     motifArr = map(lambda x: map(lambda y: y/(1.0 * n), x), motifArr)
#     return motifArr

# def savePSSMMode(motifs, n, mode, filename, likelihood):
#     motifArrays = map(lambda x: getPSSM(x), motifs)
    
#     f = open(filename, "w")
#     f.write("\nNumber of sequences: " + str(n) + "\nNumber of modes: " + str(mode) + "\nLikelihood: " + str(likelihood) + "\n\n")
#     for i in range(mode):
#         f.write("\nMODE: " + str(i) + " (" + str(len(motifs[i])) + " sequences)\n\n")
#         if motifArrays[i] is None: continue
#         f.write("PO\tA\tC\tG\tT\n");
#         for j in range(len(motifArrays[i])):
#             f.write(str(j).zfill(2) + "\t" + str(motifArrays[i][j][0]) + "\t" + str(motifArrays[i][j][1]) + "\t" + str(motifArrays[i][j][2]) + "\t" + str(motifArrays[i][j][3]) + "\n")
#         f.write("\n\n")

def savePSSM(motifs, n, minMode, maxMode, dirname, likelihoods):
    for i in range(minMode, maxMode + 1):
        savePSSMMode(motifs[i - minMode], n, i, dirname + "/mode_" + str(i) + "/pssm.txt", likelihoods[i - minMode])

def saveLogosMode(dataFile, trainOut, mode, dirname):
    # data = []
    # # j = 0
    # # n = 0
    # tmpData = ""
    # with open(dataFile) as infile:
    #     for line in infile:
    #         if line[0] == '>':
    #             if tmpData != "": data.append(tmpData)
    #             tmpData = ""
    #         else:
    #             tmpData = tmpData + line.strip()
    # if tmpData != "": data.append(tmpData)
    #         # if(j%2 == 1): 
    #         #     data.append(line[:-1])
    #         #     # n = n + 1
    #         # j = j + 1

    data = readData(dataFile)
    n = len(trainOut['labels'])
    start = []
    label = []
    width = []
    for i in range(n):
        label.append(trainOut['labels'][i])
        start.append(trainOut['startPos'][i])
    for i in range(mode):
        width.append(trainOut['motifWidth'][i])
    for i in range(n):
        ldi = len(data[i])
        if start[i] == -1: label[i] = -1
        elif start[i] < ldi:
            tmp = data[i][start[i]: (start[i] + width[label[i]])]
            data[i] = tmp
        else:
            tmpL = map(lambda x: getN(x), (data[i][::-1])[(start[i] - ldi - 1):(start[i] + width[label[i]] - ldi - 1)])
            tmpL = ''.join(tmpL)
            data[i] = tmpL

    dataList = []
    for i in range(mode): dataList.append([])
    for i in range(n):
        if label[i] == -1: continue
        dataList[label[i]].append(data[i])
    for i in range(mode): 
        createLogo(dataList[i], dirname + "/logo_" + str(i) + ".png")
        createLogo(map(lambda x: reverseSequence(x), dataList[i]), dirname + "/logo_rev_" + str(i) + ".png")
    return dataList

def saveInfoFileMode(dataFile, trainOut, mode, filename):
    # data = []
    seqName = []
    # j = 0
    # n = 0

    with open(dataFile) as infile:
        for line in infile:
            if line[0] == '>':
                tmp = re.search(r'chr[0-9a-zA-Z]+:[0-9]+-[0-9]+', line[1:])
                if tmp: seqName.append((tmp.group(0)).strip())
                else: seqName.append((line[1:]).strip())
                
    # with open(dataFile) as infile:
    #     for line in infile:
    #         if(j%2 == 1):
    #             data.append(line[:-1])
    #             # n = n + 1
    #         else:
    #             tmp = re.search(r'chr[0-9a-zA-Z]+:[0-9]+-[0-9]+', line[1:])
    #             if tmp: seqName.append((tmp.group(0)).strip())
    #             else: seqName.append((line[1:]).strip())
    #         j = j + 1

    data = readData(dataFile)
    n = len(trainOut['labels'])
    start = []
    label = []
    width = []
    strand = []

    # print "N", n
    # print "Labels length", len(trainOut['labels'])

    for i in range(n):
        label.append(trainOut['labels'][i])
        start.append(trainOut['startPos'][i])
    for i in range(mode):
        width.append(trainOut['motifWidth'][i])
    strand = []
    checkFormat = re.search(r'chr[0-9a-zA-Z]+:[0-9]+-[0-9]+', seqName[0]) is not None
    for i in range(n):
        ldi = len(data[i])
        if start[i] == -1: label[i] = -1
        elif start[i] < ldi:
            tmp = data[i][start[i]: (start[i] + width[label[i]])]
            data[i] = tmp
            strand.append("+")
            if checkFormat: start[i] = int(str.split((str.split(seqName[i], ":")[1]), "-")[0]) + start[i]
        else:
            tmpL = map(lambda x: getN(x), (data[i][::-1])[(start[i] - ldi - 1):(start[i] + width[label[i]] - ldi - 1)])
            tmpL = ''.join(tmpL)
            data[i] = tmpL
            start[i] = start[i] + width[label[i]] - ldi - 1
            strand.append("-")
            if checkFormat: start[i] = int(str.split((str.split(seqName[i], ":")[1]), "-")[0]) + ldi - start[i] 

    f = open(filename, "w")
    f.write("#sequenceName\tmodeNumber\tpositionInSequence\tstrand\tsite\n")

    j = 0
    for i in range(n):
        if label[i] == -1: f.write(seqName[i] + "\t-1\t-1\t-1\t-1\n")
        else:
            f.write(seqName[i] + "\t" + str(label[i]) + "\t" + str(start[i]) + "\t" + strand[j] + "\t" + data[i] + "\n")
            j = j + 1
    f.close()

def saveInfoFiles(d, trainOut):
    for i in range(d['-maxMode'] - d['-minMode'] + 1):
        saveInfoFileMode(d['-f'], trainOut[i], i + d['-minMode'], d['-o'][1] + "/mode_" + str(i + d['-minMode']) + "/info.txt")

def saveLogos(d, trainOut):
    motifs = []
    for i in range(d['-maxMode'] - d['-minMode'] + 1):
        tmp = saveLogosMode(d['-f'], trainOut[i], i + d['-minMode'], d['-o'][1] + "/mode_" + str(i + d['-minMode']))
        motifs.append(tmp)
    return motifs

def removeMultipleLikes(startMode, endMode, trials, dirname, lcount):
    for i in range(startMode, endMode + 1):
        os.system("cp " + dirname + "/mode_" + str(i) + "/likes_" + str(trials[i - startMode]) + ".txt" + " " + dirname + "/mode_" + str(i) + "/likes.txt")
        for j in range(lcount):
            os.system("rm -f " + dirname + "/mode_" + str(i) + "/likes_" + str(j + 1) + ".txt")

def removeMultipleLikesMode(mode, seed, dirname, lcount):
    #    for i in range(startMode, endMode + 1):
    os.system("cp " + dirname + "/mode_" + str(mode) + "/likes_" + str(seed) + ".txt" + " " + dirname + "/mode_" + str(mode) + "/likes.txt")
    for j in range(lcount):
        os.system("rm -f " + dirname + "/mode_" + str(mode) + "/likes_" + str(j + 1) + ".txt")

def saveBinaryOut(trainOut, filename):
    l = [(x['likelihood'], x['motifWidth']) for x in trainOut]
    o = pickle.dump(l, open(filename, "wb"))

def saveModeDetails(d, trainOut, seed, mode):
    os.system("cp -r " + d['-o'][1] + "/mode_" + str(mode) + "/Trial_" + str(seed) + "/* " + d['-o'][1] + "/mode_" + str(mode) + "/")
    # if d['-v'] != 0:
    #     removeMultipleLikesMode(mode, seed, d['-o'][1], d['-lcount'])
    #     plotFigures.plotLikelihoodMode(d, mode)
    # saveInfoFileMode(d['-f'], trainOut, mode, d['-o'][1] + "/mode_" + str(mode) + "/info.txt")
    # motifs = saveLogosMode(d['-f'], trainOut, mode, d['-o'][1] + "/mode_" + str(mode))
    # savePSSMMode(motifs, len(trainOut['labels']), mode, d['-o'][1] + "/mode_" + str(mode) + "/pssm.txt", trainOut['likelihood'])

def saveModeTrialDetails(d, trainOut, trial, mode):
    # if d['-v'] != 0:
    #     removeMultipleLikesMode(mode, seed, d['-o'][1], d['-lcount'])
    #     plotFigures.plotLikelihoodMode(d, mode)
    if d['-v'] != 0: plotFigures.plotSingleFile(d, d['-o'][1] + "/mode_" + str(mode) + "/Trial_" + str(trial))
    saveInfoFileMode(d['-f'], trainOut, mode, d['-o'][1] + "/mode_" + str(mode) + "/Trial_" + str(trial) + "/info.txt")
    motifs = saveLogosMode(d['-f'], trainOut, mode, d['-o'][1] + "/mode_" + str(mode) + "/Trial_" + str(trial))
    savePSSMMode(motifs, len(trainOut['labels']), mode, d['-o'][1] + "/mode_" + str(mode) + "/Trial_" + str(trial) + "/pssm.txt", trainOut['likelihood'])
    
def saveDetails(d, to):
    trainOut, modelSeeds = zip(*to)
    # likelihoods = [x['likelihood'] for x in trainOut]
    # if d['-v'] != 0: 
    #     removeMultipleLikes(d['-minMode'], d['-maxMode'], modelSeeds, d['-o'][1], d['-lcount'])
    #     plotFigures.plotLikelihood(d)
    # saveInfoFiles(d, trainOut)
    # motifs = saveLogos(d, trainOut)
    # savePSSM(motifs, len(trainOut[0]['labels']), d['-minMode'], d['-maxMode'], d['-o'][1], likelihoods)
    saveBinaryOut(trainOut, d['-o'][1] + "/models.bin.p")
