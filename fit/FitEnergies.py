#!/usr/bin/python3.6
from headerSimAnn import *

#input
BASE_NAME = "trn" ##name of split files without the final extensions
LATTICE_FILE_NAME = "lattice.tsam"
NRG_OUTPUT_NAME = "output.csv" ##name of result from writeOutput.sh
#output
ENV_OUTPUT_NAME = "fitEnergies.tsam"


#Read the energy CSV file, expects element 1 to be the ID number and element 2 to be the energy
#of that index number
def SetB(datLoc, nRows):
    res = [None for i in range(0, nRows)]
    with open(datLoc, 'r') as infile:
        for n, lin in enumerate(infile):
            if(n == 0):
                continue
            line = lin.split(',')
            res[int(line[1])] = float(line[2])
        infile.close()
    return res

A, sites, map = ReadLatticeFile(LATTICE_FILE_NAME) ##note the A is the 3x3 lattice vector matrix, NOT
envs = ReadEnvFile(infileLoc=BASE_NAME + ".env")   ##the usual 'A' in the matrix equation that we're
dec = ReadDecompFile(BASE_NAME + ".dec")           ##solving (Ax ~ b).
occs = ReadOccFile(BASE_NAME + ".occ")
b = SetB(datLoc=NRG_OUTPUT_NAME, nRows=len(dec))


##Do the solving.  Make sure the result, x, is a python list
##'dec' is the coefficient matrix and b contains the energies of the cells
"""Here, you should do your own fitting routine if you don't have scipy"""
"""I like scipy's sparse lsqr because it allows for a physically relavent initial guess"""
"""(and also, we do have a sparse matrix...)"""
"""If you use numpy's least squares routine, for example, your "energies" could be well over 1E6!"""
from numpy import array
from scipy.sparse.linalg import lsqr

x0 = [sum(b)/len(b)/len(occs[0]) for i in range(0, len(dec[0]))]  ##initial guess: avg embedding nrg / site
x = list( lsqr(A=array(dec), b=array(b), x0=array(x0))[0] )
"""No need to change anything past this point"""


#Write results to output file
for i in range(0, len(x)):
    envs[i].SetEnergy(x[i])
WriteEnvFile(ENV_OUTPUT_NAME, envs)
