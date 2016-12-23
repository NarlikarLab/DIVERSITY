import os

def plotSingleFile(d, dirname):
    f1 = dirname + "/likes.txt"
    f2 = dirname + "/likes.png"
    os.system("gnuplot -e 'filename=\"" + f1 + "\"; var=\"" + f2 + "\"' " + d['-v'])
    os.system("rm " + f1)

def plotLikelihood(d):
    for i in range(d['-minMode'], d['-maxMode'] + 1):
        f1 = d['-o'][1] + "/mode_" + str(i) + "/likes.txt"
        f2 = d['-o'][1] + "/mode_" + str(i) + "/likes.png"
        os.system("gnuplot -e 'filename=\"" + f1 + "\"; var=\"" + f2 + "\"' " + d['-v'])
        os.system("rm " + f1)


def plotLikelihoodMode(d, mode):
    f1 = d['-o'][1] + "/mode_" + str(mode) + "/likes.txt"
    f2 = d['-o'][1] + "/mode_" + str(mode) + "/likes.png"
    os.system("gnuplot -e 'filename=\"" + f1 + "\"; var=\"" + f2 + "\"' " + d['-v'])
    os.system("rm " + f1)
