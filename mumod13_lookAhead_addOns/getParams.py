import re
import sys
import os

def printHelp():
    print "Usage:\n./mumod [options]"
    print "\t-f fasta file(compulsory)"
    print "\t-o output prefix"
    print "\t-a pseudo count. Default 1"
    print "\t-maskReps mask repeats in sequences. Default 1(mask)"
    print "\t-r include reverse strand while training. Default 1(include)"
    print "\t-zoops zero or one occurence per sequence. 0 means all sequences must have a motif. 1 means all sequences may not have a motif. Any value between 0 and 1 implies the probability of a sequence not having a motif. Default 0"
    print "\t-minWidth minimum motif width. Default 6"
    print "\t-initialWidth starting width of motifs. Default 8"
    print "\t-minMode minimum number of modes. Default 1"
    print "\t-maxMode maximum number of modes. Default 10"
    print "\t-lcount number of models to be learned while training. Best model is considered. Default 5"
    print "\t-proc maximum number of processing units to use for computation. Default is the number of processing units the system has"
    print "\t-v 0 or 1. 1 to save likelihood plots of learned models. Default 0"
    exit(2)


def validInt(s, opt):    # Check if string is a valid positive integer
    s1 = re.search('^[1-9][0-9]*\Z', s)
    if s1 is None: 
        print "ERROR: Invalid option: " + opt + " " + s + ". Must be positive integer"
        printHelp()
    s = s1
    return int(s.group(0))

def validFD(s, opt):    # Check is string is a valid filename
    s1 = re.search('^.*\Z', s)
    if s1 is None: 
        print "ERROR: Invalid filename " + s
        printHelp()
    s = s1
    return s.group(0)

def validR(s, opt):    # Check if string is 0 or 1
    s1 = re.search('^(0|1)\Z', s)
    if s1 is None: 
        print "ERROR: Invalid option: " + opt + " " + s + ". Must be 0 or 1"
        printHelp()
    s = s1
    s = int(s.group(0))
    return s

def validNum(s, opt):    # Check is string is a valid positive number
    s1 = re.search('^[0-9]+?(\.[0-9]+)?\Z', s)
    if s1 is None: 
        print "ERROR: Invalid option: " + opt + " " + s + ". Must be a positive real number"
        printHelp()
    s = s1
    s = float(s.group(0))
    if s == 0: 
        printHelp()
    return s

def validNum1(s, opt):    # Check is string is a valid non negative number less than or equal to 1
    s1 = re.search('^[0-9]+?(\.[0-9]+)?\Z', s)
    if s1 is None: 
        print "ERROR: Invalid option: " +opt + " " + s + ". Must be a real number between 0 and 1"
        printHelp()
    s = s1
    s = float(s.group(0))
    if s > 1: printHelp()
    return s

def validFile(s):    # Check if string is a valid file
    if s == "": printHelp(1)
    if not os.path.isfile(s):
        print("Could not open file: " + s)
        exit(2)
    f = open(s, "r")
    lines = f.readlines()
    f.close
    return len(lines), lines[0]

def validDir(s):    # Check if string is a valid directory
    defVal = 0
    if s == "": 
        s = "mumodOut"
        defVal = 1
    elif s[-1] == "/":
        s = s[:-1]
    old = s
    d = os.path.dirname(s)
    if d != "" and not os.path.isdir(d):
        print("Invalid directory: " + d)
        exit(2)
    if os.path.isdir(s):
        b = 1
        while True:
            if os.path.isdir(s + "_" + str(int(b))): b = b + 1
            else:
                new = s + "_" + str(int(b))
                break
    else: new = s
    return [old, new, defVal]

def getValues():
    d = {'-f': '', '-r': 1, '-a': 1, '-lambda': 0, '-zoops': 0, '-minWidth': 6, '-o': '', '-minMode': 1, '-maxMode': 10, '-initialWidth': 8, '-lcount': 5, '-proc': 0, '-maskReps': 1, '-v': 0}
    dF = {'-f': validFD, '-r': validR, '-a': validNum, '-zoops': validNum1, '-minWidth': validInt, '-o': validFD, '-minMode': validInt, '-maxMode': validInt, '-initialWidth': validInt, '-lcount': validInt, '-proc': validInt, '-maskReps': validR, '-v': validR}

    lst = []
    sysargv = (sys.argv)[1:]
    if sysargv[:-1] == []: printHelp()
    for i in range(len(sysargv)):
        if i % 2 == 0: tmp = sysargv[i]
        else: lst = lst + [tmp + " " + sysargv[i]]
    for i in lst:
        s = re.search('^\-[a-zA-Z]* ', i)
        if s is None: 
            print "ERROR: Invalid option", i.split(" ")[0]
            printHelp()
        s = s.group(0).strip()
        if s in d: d[s] = dF[s](i.split()[1], s)
        else:
            print "ERROR: Invalid option", s
            printHelp()
    if(d['-maxMode'] < d['-minMode']):
        print "ERROR: -maxMode cannot be less than -minMode"
        exit(2)
    validFile(d['-f'])
    d['-o'] = validDir(d['-o'])
    try:
        os.mkdir(d['-o'][1])
    except:
        print "ERROR: Cannot create directory", d['-o'][1]
        exit(2)

    for i in range(d['-minMode'], (d['-maxMode'] + 1)):
        try:
            os.mkdir(d['-o'][1] + "/mode_" + str(i))
            for j in range(1, (d['-lcount'] + 1)):
                try:
                    os.mkdir(d['-o'][1] + "/mode_" + str(i) + "/Trial_" + str(j))
                except:
                    "ERROR: Cannot create directory" + d['-o'][1] + "/mode_" + str(i) + "/Trial_" + str(j)
        except:
            print "ERROR: Cannot create directory" + d['-o'][1] + "/mode_" + str(i)
    if d['-v'] != 0: d['-v'] = sysargv[-1] + "makeLikelihoodPlot.pg"
    return d