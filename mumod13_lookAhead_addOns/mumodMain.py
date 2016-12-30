import multiprocessing as mp
from cstructures import *
from getParams import *
import evaluate as ev
import saveFiles as sf

# print execution information
def printDetails(d):
    print "\n\nProcessing details:\n"
    if d['-o'][2] != 0:
        print "Output directory was not provided. Output would be saved in default directory", d['-o'][1]
    else:
        if d['-o'][0] != d['-o'][1]:
            print "Directory", d['-o'][0], "was already present"
        print "Output saved in", d['-o'][1]
    print "Output file consisting of motif positions and strand information would be saved as " + d['-o'][1] + "/info.txt for each mode"
    print "Sequence logos for all motifs of every mode would be saved in each mode directory"
    print "PSSMs for all motifs of every mode would be saved as " + d['-o'][1] + "/pssm.txt"
    print "\n\n"

# function that calls the training function and then calls another function to save details
def getModel(d):
    dirname = d['-o'][1]
    n = mp.cpu_count()
    if d['-proc'] < n and d['-proc'] > 0: n = d['-proc']
    d['-proc'] = n
    printDetails(d)
    to = ev.learn(d)
    print "\n\nSaving details..."
    sf.saveDetails(d, to)
    print "\nGoodbye!"

if __name__ == '__main__':
    try:
        d = None
        d = getValues()
        getModel(d)
    except (KeyboardInterrupt, SystemExit):
        os.system('setterm -cursor on')
        if d is not None: 
            print "\nExiting...\n"
        exit(2)
    exit(0)
