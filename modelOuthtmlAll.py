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

# Create HTML file for motifs of all models

import sys
import os
from config import *

def createHTMLAll(infastafile, opdir, minMode, maxMode, bestModel):
    f = file(opdir+'/' + htmlFileAll,'w')
    f.write("<!DOCTYPE html>\n<html>\n<body>\n<h1>MODEL</h1>\n")
    f.write("<h3> Input File: "+ infastafile + " </h3>\n")
    f.write("<h3> Number of models learned: "+str(maxMode-minMode+1)+" </h3>\n")
    f.write("<h3> Number of modes in best model: "+str(bestModel)+" </h3>\n")
    logodir = modeDir.format(bestModel)

    for j in range(minMode, maxMode + 1):
        s = "(Best model)" if j == bestModel else ""
        f.write("<h4><a href=\"" + htmlFileMode.format(j) + "\">Model " + str(j) + s + "</a></h4>\n")
        logodir = modeDir.format(j)
        f.write("<h5>Forward strand motifs</h5>\n")
        f.write("<div id=\"logo_images\">\n")
        for i in range(j):
            if os.path.isfile(opdir + "/" + logodir + "/" + motifLogoName.format(i)):
                f.write("<figure>\n")
                f.write("<img height=\"50px\" width=\"200px\" src=\"" + logodir + "/" + motifLogoName.format(i) + "\">\n")
                f.write("<figcaption><i> " + "mode_" + str(i + 1) + "</i></figcaption>\n")
		f.write("</figure>\n")
            else:
                continue
        f.write("</div>\n")
        f.write("<br>\n")
        f.write("<h5>Reverse strand motifs</h5>\n")
        f.write("<div id=\"logo_images\">\n")
        for i in range(j):
            if os.path.isfile(opdir + "/" + logodir + "/" + motifLogoReverseName.format(i)):
                f.write("<figure>\n")
                f.write("<img height=\"50px\" width=\"200px\" src=\"" + logodir + "/" + motifLogoReverseName.format(i) + "\">\n")
                f.write("<figcaption><i> " + "mode_" + str(i + 1) + "</i></figcaption>\n")
		f.write("</figure>\n")
            else:
                continue
        f.write("</div>\n")
        f.write("<br>\n")

    f.write("<style>\n\t div{padding:0;}\n\t figure{display: inline-block;}\n\t img{height:auto;}\n\t")
    f.write("figure figcaption{border:35px black; text-align:center;}\n\t")
    f.write("th,td{text-align:center;border-spacing:10px;border:1px solid black;height:50px;width:100px;}\n</style>\n")
    
    f.write("<p><i>NOTE: Best model is chosen using lambda = 5</i></p>\n")
    f.write("</body>\n</html>\n")
    f.close()

if __name__ == '__main__':
    infastafile = sys.argv[1]
    opdir = sys.argv[2]
    minMode = int(sys.argv[3])
    maxMode = int(sys.argv[4])
    createHTMLAll(infastafile, opdir, minMode, maxMode, 5)
