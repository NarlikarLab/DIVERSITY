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

def chooseBestModel(dirname, lmbda):
    with open(dirname + "/" + modelBinaryFile, "rb") as f:
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

v = chooseBestModel(d, 5)
print "\nBest Model:", v

v = chooseBestModel(d, 10)
print "\nBest Model:", v

#print "\nBest Model:", chooseBestModel(d, 10)

# for i in range(1, 50):
#     print i, chooseBestModel(d, i)
