import pickle
import os
import sys

def saveBestModel(dirname, lmbda, mode):
    os.system("cp -r " + dirname + "/mode_" + str(mode) + " " + dirname + "/bestModelLambda" + str(lmbda))

def chooseBestModel(dirname1, dirname2, lmbda):
    with open(dirname1 + "/models.bin.p", "rb") as f:
        l1 = pickle.load(f)
    likes1 = [x[0] - 3*lmbda*sum(x[1]) for x in l1]
    with open(dirname2 + "/models.bin.p", "rb") as f:
        l2 = pickle.load(f)
    likes2 = [x[0] - 3*lmbda*sum(x[1]) for x in l2]
    index1 = likes1.index(max(likes1))
    index2 = likes2.index(max(likes2))

    if likes1[index1] > likes2[index2]:
        mode = len(l1[index1][1])
        it = 1
    else:
        mode = len(l2[index2][1])
        it = 2
    # print "iter 1:", likes1[index1], "iter 2:", likes2[index2]
    return mode, it
    # mode = len(l[index][1])
    # # saveBestModel(dirname, lmbda, mode)
    # return len(l[index][1])

# chooseBestModel((sys.argv)[1], 5)

d1 = (sys.argv)[1]
d2 = (sys.argv)[2]
for i in range(1, 10):
    print i, chooseBestModel(d1, d2, i)
