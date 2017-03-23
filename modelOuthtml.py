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

# Create HTML file

import sys
import os
from config import *

def createHTML(infastafile, opdir, minMode, maxMode, bestModel):
    f = file(opdir+'/' + htmlFile,'w')
    f.write("<!DOCTYPE html>\n<html>\n<body>\n<h1>MODEL</h1>\n")
    f.write("<h3> Input File: "+ infastafile + " </h3>\n")
    f.write("<h3> Number of models learned: "+str(maxMode-minMode+1)+" </h3>\n")
    f.write("<h3> Number of modes in best model: "+str(bestModel)+" </h3>\n")
    logodir = modeDir.format(bestModel)

    for i in range(bestModel):
        f.write("<h4> Mode " + str(i + 1) + "</h4>\n")
        f.write("<img src=\"" + logodir + "/" + motifLogoName.format(i) + "\" style=\"border:thin solid black\">\n")

    f.write("<p><i>NOTE: Best model is chosen using lambda = 5</i></p>\n")
    f.write("</body>\n</html>\n")
    f.close()

if __name__ == '__main__':
    infastafile = sys.argv[1]
    opdir = sys.argv[2]
    minMode = int(sys.argv[3])
    maxMode = int(sys.argv[4])
    createHTML(infastafile, opdir, minMode, maxMode, 5)
