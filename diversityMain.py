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

# Main function executing DIVERSITY

import multiprocessing as mp
from cstructures import *
from getParams import *
import evaluate as ev
import saveFiles as sf

# save settings in file

def saveSettings(d):
    f = open(d['-o'][1] + "/" + settingsFile, "w")
    f.write("\nInput fasta file: " + os.path.abspath(d['-f']))
    f.write("\n\nOutput prefix: " + os.path.abspath(d['-o'][1]))
    f.write("\n\nMinimum number of modes: " + str(d['-minMode']))
    f.write("\n\nMaximum number of modes: " + str(d['-maxMode']))
    f.write("\n\nNumber of models to be learnt per mode: " + str(d['-lcount']))
    f.write("\n\nMinimum motif width: " + str(d['-minWidth']))
    f.write("\n\nMaximum motif width: " + str(d['-maxWidth']))
    f.write("\n\nInitial motif width: " + str(d['-initialWidth']))
    f.write("\n\nMask repeats: ")
    if (d['-maskReps'] == 1): f.write("Yes")
    else: f.write("No")
    f.write("\n\nInclude reverse strand: ")
    if (d['-r'] == 1): f.write("Yes")
    else: f.write("No")
    f.write("\n\nZOOPS: " + str(d['-zoops']))
    f.write("\n\nMaximum number of cores used: " + str(d['-proc']))
    f.write("\n")
    f.close()
    
# print execution information
def printDetails(d):
    print "\n\nProcessing details:\n"
    if d['-o'][2] != 0:
        print "Output directory was not provided. Output would be saved in default directory", d['-o'][1]
    else:
        if d['-o'][0] != d['-o'][1]:
            print "Directory", d['-o'][0], "was already present"
        print "Output saved in", d['-o'][1]
    print "Output file consisting of motif positions and strand information would be saved as " + d['-o'][1] + "/" + infoFile + " for each mode"
    print "Sequence logos for all motifs of every mode would be saved in each mode directory"
    print "PSSMs for all motifs of every mode would be saved as " + d['-o'][1] + "/" + pssmFile
    print("Execution settings will be saved as " + d['-o'][1] + "/" + settingsFile)
    print "\n\n"

# function that calls the training function and then calls another function to save details
def getModel(d):
    dirname = d['-o'][1]
    n = mp.cpu_count()
    if d['-proc'] < n and d['-proc'] > 0: n = d['-proc']
    d['-proc'] = n
    printDetails(d)
    saveSettings(d)
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
