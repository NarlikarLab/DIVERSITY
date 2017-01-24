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

import pickle
import os
import sys
from config import *

def saveBestModel(dirname, lmbda, mode):
    os.system("cp -r " + dirname + "/" + modeDir.format(str(mode)) + " " + dirname + "/bestModelLambda" + str(lmbda))

def chooseBestModel(dirname1, dirname2, lmbda):
    with open(dirname1 + "/" + modelBinaryFile, "rb") as f:
        l1 = pickle.load(f)
    likes1 = [x[0] - 3*lmbda*sum(x[1]) for x in l1]
    with open(dirname2 + "/" + modelBinaryFile, "rb") as f:
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
