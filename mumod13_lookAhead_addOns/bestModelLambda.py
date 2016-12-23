import pickle
import os
import sys

def saveBestModel(dirname, lmbda, mode):
    os.system("cp -r " + dirname + "/mode_" + str(mode) + " " + dirname + "/bestModelLambda" + str(lmbda))

def chooseBestModel(dirname, lmbda):
    with open(dirname + "/models.bin.p", "rb") as f:
        l = pickle.load(f)
#    print [len(x[1]) for x in l]
    likes = [x[0] - 3*lmbda*sum(x[1]) for x in l]
#    print likes
    index = likes.index(max(likes))
    print "\n\nLambda:", lmbda
    print "#modes likelihood penalized"
    for i in range(len(likes)):
        print len(l[i][1]), l[i][0], likes[i]
    mode = len(l[index][1])
    # saveBestModel(dirname, lmbda, mode)
    return len(l[index][1])

# chooseBestModel((sys.argv)[1], 5)

d = (sys.argv)[1]

v = chooseBestModel(d, 3)
print "\nBest Model:", v

v = chooseBestModel(d, 5)
print "\nBest Model:", v

v = chooseBestModel(d, 10)
print "\nBest Model:", v

#print "\nBest Model:", chooseBestModel(d, 10)

# for i in range(1, 50):
#     print i, chooseBestModel(d, i)
