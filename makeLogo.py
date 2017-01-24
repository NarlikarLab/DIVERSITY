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


import weblogoMod.weblogolib as wl
import sys
from config import *

# function to create a logo using modified Weblogo3.3 from given sequences
def createLogo(sequences, filename):
    seqs = wl.read_seq_data(sequences)
    data = wl.LogoData.from_seqs(seqs)
    options = wl.LogoOptions()
    options.title = ""
    options.size = "large"
    options.color_scheme = wl.colorscheme.monochrome
    formt = wl.LogoFormat(data, options)
    fout = open(filename, "w")
    wl.png_formatter(data, formt, fout)
    fout.close()

# get nucleotide for reverse strand
def getN(c):
    if c == 'N': return c
    elif c == 'A': return 'T'
    elif c == 'C': return 'G'
    elif c == 'G': return 'C'
    else: return 'A'

def getLogos():
    dataFile = (sys.argv)[1]
    labelFile = (sys.argv)[2]
    startFile = (sys.argv)[3]
    widthFile = (sys.argv)[4]

    data = []
    j = 0
    n = 0
    with open(dataFile) as infile:
        for line in infile:
            if(j%2 == 1): 
                data.append(line[:-1])
                n = n + 1
            j = j + 1

    label = []
    with open(labelFile) as infile:
        for line in infile:
            label.append(int(line))

    start = []
    with open(startFile) as infile:
        for line in infile:
            start.append(int(line))

    width = []
    with open(widthFile) as infile:
        for line in infile:
            width.append(int(line))

    for i in range(n):
        ldi = len(data[i])
        if start[i] < ldi:
            tmp = data[i][start[i]: (start[i] + width[label[i]])]
            data[i] = tmp
        else:
            tmpL = map(lambda x: getN(x), (data[i][::-1])[(start[i] - ldi - 1):(start[i] + width[label[i]] - ldi - 1)])
            tmpL = ''.join(tmpL)
            data[i] = tmpL

    modes = max(label) + 1
    dataList = []
    for i in range(modes): dataList.append([])
    for i in range(n):
        dataList[label[i]].append(data[i])

    for i in range(modes): 
        print "Width", width[i]
        createLogo(dataList[i], modeDir.format(str(i)) + ".png")

getLogos()
