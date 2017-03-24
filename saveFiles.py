##################### DIVERSITY #####################

#    DIVERSITY is a tool to explore multiple ways of protein-DNA
#    binding in the genome. More information can be found in the README file.
#    Copyright (C) 2015  Sneha Mitra, Anushua Biswas and Leelavati Narlikar

#    DIVERSITY is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.

#    DIVERSITY is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.

#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.

######################################################

# Save DIVERSITY output files

import weblogoMod.weblogolib as wl
import numpy as np
import pickle
import os
import re
import plotFigures
from config import *
from bestModelLambda import chooseBestModel
from modelOuthtml import createHTML
from modelOuthtmlAll import createHTMLAll

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

# function to return nucleotide in reverse strand
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

# map A, C, G, T to numbers
def nToNum(c):
    if c == 'A' or c == 'a': return 0
    elif c == 'C' or c == 'c': return 1
    elif c == 'G' or c == 'g': return 2
    else: return 3

# get reverse strand sequence    
def reverseSequence(sequence):
    return ''.join(map(lambda x: getN(x), list(sequence))[::-1])

# create logo for given set of sequences using a modifed version of Weblogo3.3
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

# generate PSSM for given motif
def getPSSM(motif):
    n = len(motif)
    if n == 0: return
    mWidth = len(motif[0])
    motifArr = np.array(map(lambda x: map(lambda y: nToNum(y), list(x)), motif)).transpose()
    motifArr = map(lambda x: np.bincount(x, minlength = 4), motifArr)
    motifArr = map(lambda x: map(lambda y: y/(1.0 * n), x), motifArr)
    motifArr = np.array(motifArr).transpose()
    return motifArr

# save PSSMs for a given model
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

# save PSSMs for all models
def savePSSM(motifs, n, minMode, maxMode, dirname, likelihoods):
    for i in range(minMode, maxMode + 1):
        savePSSMMode(motifs[i - minMode], n, i, dirname + "/" + modeDir.format(str(i)) + "/" + pssmFile, likelihoods[i - minMode])

# save logos for a given model
def saveLogosMode(dataFile, trainOut, mode, dirname):
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
        createLogo(dataList[i], dirname + "/" + motifLogoName.format(str(i)))
        createLogo(map(lambda x: reverseSequence(x), dataList[i]), dirname + "/" + motifLogoReverseName.format(str(i)))
    return dataList

# save model information for a given model
def saveInfoFileMode(dataFile, trainOut, mode, likesInfoFile, filename):
    seqName = []
    with open(dataFile) as infile:
        for line in infile:
            if line[0] == '>':
                tmp = re.search(r'chr[0-9a-zA-Z]+:[0-9]+-[0-9]+', line[1:])
                if tmp: seqName.append((tmp.group(0)).strip())
                else: seqName.append((line[1:]).strip())

    data = readData(dataFile)
    n = len(trainOut['labels'])
    start = []
    label = []
    width = []
    strand = []

    pMode = []
    likes = []
    with open(likesInfoFile) as infile:
        for line in infile:
            l = line.strip().split()
            l = map(lambda x: float(x), l)
            if not pMode: pMode = l
            else: likes.append(l)
    
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
    f.write("#sequenceName\tmodeNumber\tpositionInSequence\tstrand\tsite\tScore\n")

    j = 0
    for i in range(n):
        l = str(sum(map(lambda x: likes[i][x]*pMode[x], range(mode))))
        if label[i] == -1: f.write(seqName[i] + "\t-1\t-1\t-1\t-1\t0\t" + l + "\n")
        else:
            f.write(seqName[i] + "\t" + str(label[i]) + "\t" + str(start[i]) + "\t" + strand[j] + "\t" + data[i] + "\t" + l + "\n")
            j = j + 1
    f.close()
    os.system("rm -f " + likesInfoFile)

# save model information for all learned models
def saveInfoFiles(d, trainOut):
    for i in range(d['-maxMode'] - d['-minMode'] + 1):
        saveInfoFileMode(d['-f'], trainOut[i], i + d['-minMode'], d['-o'][1] + "/" + modeDir.format(str(i + d['-minMode'])) + "/" + infoFile)

# save logos for all learned models
def saveLogos(d, trainOut):
    motifs = []
    for i in range(d['-maxMode'] - d['-minMode'] + 1):
        tmp = saveLogosMode(d['-f'], trainOut[i], i + d['-minMode'], d['-o'][1] + "/" + modeDir.format(str(i + d['-minMode'])))
        motifs.append(tmp)
    return motifs


# remove redundant files
def removeMultipleLikes(startMode, endMode, trials, dirname, lcount):
    for i in range(startMode, endMode + 1):
        os.system("cp " + dirname + "/" + modeDir.format(str(i)) + "/" + likelihoodTrialFile.format(str(trials[i - startMode])) + " " + dirname + "/" + modeDir.format(str(i)) + "/" + likelihoodFile)
        for j in range(lcount):
            os.system("rm -f " + dirname + "/" + modeDir.format(str(i)) + "/" + likelihoodTrialFile.format(str(j + 1)))

# remove redundant files
def removeMultipleLikesMode(mode, seed, dirname, lcount):
    os.system("cp " + dirname + "/" + modeDir.format(str(mode)) + "/" + likelihoodTrialFile.format(str(seed)) + " " + dirname + "/" + modeDir.format(str(mode)) + "/" + likelihoodFile)
    for j in range(lcount):
        os.system("rm -f " + dirname + "/" + modeDir.format(str(mode)) + "/" + likelihoodTrialFile.format(str(j + 1)))

# save model in binary format in python pickle form
def saveBinaryOut(trainOut, filename):
    l = [(x['likelihood'], x['motifWidth']) for x in trainOut]
    o = pickle.dump(l, open(filename, "wb"))

# save model information
def saveModeDetails(d, trainOut, seed, mode):
    os.system("cp -r " + d['-o'][1] + "/" + modeDir.format(str(mode)) + "/" + trialDir.format(str(seed)) + "/* " + d['-o'][1] + "/" + modeDir.format(str(mode)) + "/")

# save model information for all trials
def saveModeTrialDetails(d, trainOut, trial, mode):
    if d['-v'] != 0: plotFigures.plotSingleFile(d, d['-o'][1] + "/" + modeDir.format(str(mode)) + "/" + trialDir.format(str(trial)))
    saveInfoFileMode(d['-f'], trainOut, mode, d['-o'][1] + "/" + modeDir.format(str(mode)) + "/" + trialDir.format(str(trial)) + "/" + temporaryInfoFile, d['-o'][1] + "/" + modeDir.format(str(mode)) + "/" + trialDir.format(str(trial)) + "/" + infoFile)
    motifs = saveLogosMode(d['-f'], trainOut, mode, d['-o'][1] + "/" + modeDir.format(str(mode)) + "/" + trialDir.format(str(trial)))
    savePSSMMode(motifs, len(trainOut['labels']), mode, d['-o'][1] + "/" + modeDir.format(str(mode)) + "/" + trialDir.format(str(trial)) + "/" + pssmFile, trainOut['likelihood'])

# save best model in HTML format
def saveHTML(d):
    bestModel = chooseBestModel(d['-o'][1], defaultLambda)
    createHTML(d['-f'], d['-o'][1], d['-minMode'], d['-maxMode'], bestModel)
    createHTMLAll(d['-f'], d['-o'][1], d['-minMode'], d['-maxMode'], bestModel)
    
def saveDetails(d, to):
    trainOut, modelSeeds = zip(*to)
    saveBinaryOut(trainOut, d['-o'][1] + "/" + modelBinaryFile)
    saveHTML(d)
    os.system("rm -f " + d['-o'][1] + "/" + temporaryDataFile)
    os.system("rm -f " + d['-o'][1] + "/" + temporaryBinaryFile)
    
