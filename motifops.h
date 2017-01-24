///////////////////// DIVERSITY /////////////////////

//    DIVERSITY is a tool to explore multiple ways of protein-DNA
//    binding in the genome. More information can be found in the README file.
//    Copyright (C) 2015  Sneha Mitra, Anushua Biswas and Leelavati Narlikar

//    DIVERSITY is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.

//    DIVERSITY is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.

//    You should have received a copy of the GNU General Public License
//    along with this program.  If not, see <http://www.gnu.org/licenses/>.

//////////////////////////////////////////////////////

#ifndef _motifops_h
#define _motifops_h

double motifAddLeftScore(dataSet*, model*, double**, int*, int*, int);
double motifRemoveLeftScore(dataSet*, model*, double**, int*, int*, int);
double motifAddRightScore(dataSet*, model*, double**, int*, int*, int);
double motifRemoveRightScore(dataSet*, model*, double**, int*, int*, int);
void motifLeftIncrease(dataSet*, model*, double**, int*, int*, int);
void motifLeftDecrease(dataSet*, model*, double**, int*, int*, int);
void motifRightIncrease(dataSet*, model*, double**, int*, int*, int);
void motifRightDecrease(dataSet*, model*, double**, int*, int*, int);

#endif
