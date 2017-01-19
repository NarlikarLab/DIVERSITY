import os
from config import *

# plot a single likelihood file
def plotSingleFile(d, dirname):
    f1 = dirname + "/" + likelihoodFile
    f2 = dirname + "/" + likelihoodPlotFile
    os.system("gnuplot -e 'filename=\"" + f1 + "\"; var=\"" + f2 + "\"' " + d['-v'])
    os.system("rm " + f1)

# plot likelihood for all modes in different files
def plotLikelihood(d):
    for i in range(d['-minMode'], d['-maxMode'] + 1):
        f1 = d['-o'][1] + "/" + modeDir.format(str(i)) + "/" + likelihoodFile
        f2 = d['-o'][1] + "/" + modeDir.format(str(i)) + "/" + likelihoodPlotFile
        os.system("gnuplot -e 'filename=\"" + f1 + "\"; var=\"" + f2 + "\"' " + d['-v'])
        os.system("rm " + f1)


def plotLikelihoodMode(d, mode):
    f1 = d['-o'][1] + "/" + modeDir.format(str(mode)) + "/" + likelihoodFile
    f2 = d['-o'][1] + "/" + modeDir.format(str(mode)) + "/" + likelihoodPlotFile
    os.system("gnuplot -e 'filename=\"" + f1 + "\"; var=\"" + f2 + "\"' " + d['-v'])
    os.system("rm " + f1)
