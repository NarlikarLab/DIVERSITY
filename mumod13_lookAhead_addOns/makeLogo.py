import weblogoMod.weblogolib as wl
import sys

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
        createLogo(dataList[i], "mode_" + str(i) + ".png")

getLogos()
