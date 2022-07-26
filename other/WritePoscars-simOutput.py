#!/usr/bin/python3.6

# *** Meant to parse the output of the simulated annealing program ***
#Use:
#./WritePoscar.py -l <lattice file loc> -i <rsim output loc> -o <output poscar files loc>
#                 -x <<val> <val element strings to exclude from poscar>>
#                 -c <comment>
#Defaults:
#<lattice file loc> = lattice.tsam
#<rsim output loc> = out
#<output poscar file loc> = POSCAR
#<<val> <val element strings to exclude from poscar>> = 2 Va Vac
#<comment> = sa min

from sys import argv
import os
from headerSimAnn import *

#Parse CLAs
LATTICE_FILE_NAME = "lattice.tsam"
SA_CODE_FULL_OUTPUT = "out"
OUTFILE_NAME = "POSCAR"
EXCLUDE = ["Va", "Vac"]
COMMENT = "sa min"
for i in range(0, len(argv)):
    if(argv[i] == "-l"):
        LATTICE_FILE_NAME = argv[i + 1]
        continue
    if(argv[i] == "-i"):
        SA_CODE_FULL_OUTPUT = argv[i + 1]
        continue
    if(argv[i] == "-o"):
        OUTFILE_NAME = argv[i + 1]
        continue
    if(argv[i] == "-x"):
        EXCLUDE = [None for j in range(0, int(argv[i + 1]))]
        for j in range(0, int(argv[i + 1])):
            EXCLUDE[j] = argv[i + 2 + j]
        continue
    if(argv[i] == "-c"):
        COMMENT = argv[i + 1]
        continue


#Parse sim anneal output
nSites, occs, nReads = None, [], 0
with open(SA_CODE_FULL_OUTPUT, 'r') as infile:
    read = False
    for line in infile:
        if("JOB: BEGIN OPT SITE PRINTOUT" in line):
            read = True
        if(read):
            s = (infile.readline()).split()
            nReads = int(s[0])
            nSites = int(s[1])

            for i in range(0, nReads):
                occsS = (infile.readline()).split()
                occs.append([int(o) for o in occsS])
            break
    infile.close()

#Read in sites, mapping, ...
A, sites, map = ReadLatticeFile(LATTICE_FILE_NAME)

#Print POSCARs
for i in range(0, nReads):
    if(os.path.exists(str(i))):
        print("dir " + str(i) + " already exists...\n" \
              "I REFUSE TO CONTINUE WITH THIS SICK JOB ... BYE!!!")
        exit(1)
    os.mkdir(str(i))
    WritePoscar_q(str(i) + "//" + OUTFILE_NAME, A=A, sites=sites, occs=occs[i], map=map, exclude=EXCLUDE,
                  comment=COMMENT)
