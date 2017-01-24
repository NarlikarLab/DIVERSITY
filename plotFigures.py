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

# Plot likelihood values

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
