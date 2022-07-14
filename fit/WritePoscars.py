#!/usr/bin/python3.6
import os
from headerSimAnn import *

#input
BASE_NAME = "ranDecRes4" ##name of split files without the final extensions
LATTICE_FILE_NAME = "lattice-agbii4.tsam"
OFFSET = 0
SPECS_TO_EXCLUDE = ["Va", "Vac"]

A, sites, map = ReadLatticeFile(LATTICE_FILE_NAME)
envs = ReadEnvFile(infileLoc=BASE_NAME + ".env")
dec = ReadDecompFile(BASE_NAME + ".dec")
occs = ReadOccFile(BASE_NAME + ".occ")

if(not(os.path.exists("all"))):
    os.mkdir("all")

for i in range(0, len(occs)):
    #If the path already exists, this must be from a warm-start run, and we should skip it
    if(os.path.exists("all//" + str(i + OFFSET))):
        continue
    #Otherwise, go ahead and make the new directory / file
    os.mkdir("all//" + str(i + OFFSET))
    WritePoscar_q(outfileLoc="all//" + str(i + OFFSET) + "//POSCAR",
                  A=A, sites=sites, occs=occs[i], map=map, exclude=SPECS_TO_EXCLUDE,
                  comment="gen " + str(i + OFFSET))
