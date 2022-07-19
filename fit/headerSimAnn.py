#!/usr/bin/python3.6
#Header for simulated annealing stuff
DEFAULT_PREC = 11
CART_EPS = 1.0E-4

def ffmt(s, dec=DEFAULT_PREC):
    return f'{float(s):.{dec}f}'

def DistSqd_s(a, b):
    ret = 0
    for i in range(0, 3):
        ret += (b.crdsC[i] - a.crdsC[i])**(2)
    ret += (b.spec - a.spec)**(2)
    return ret

class site:
    spec = None
    crdsD = [None, None, None]
    crdsC = [None, None, None]

    def __init__(self):
        self.spec = None
        self.crdsD = [None, None, None]
        self.crdsC = [None, None, None]

    ##lineString = "2 0.25 0.65 0.89 1.0 2.0 3.0 ": 2 = species int, 3 direct coords, 3 cart crds
    def InitWithSpec(self, lineString, cCrds=True):
        split = lineString.split()
        self.spec = int(split[0])
        for i in range(0, 3):
            self.crdsD[i] = float(split[i + 1])
        if(cCrds):
            for i in range(0, 3):
                self.crdsC[i] = float(split[i + 1 + 3])

    ##linsString = "0.25 0.65 0.89": no species int: just 3 direct coords, 3 cart crds
    def InitWithoutSpec(self, lineString, cCrds=True):
        split = lineString.split()
        for i in range(0, 3):
            self.crdsD[i] = float(split[i])
        if(cCrds):
            for i in range(0, 3):
                self.crdsC[i] = float(split[i + 3])

class env:
    nSites = None
    sites = None
    nrg = None
    isRepped = None

    def __init__(self):
        self.nSites = 0
        self.sites = [] ##yes, this is necessary
        self.nrg = None
        self.isRepped = None

    def EqTo(self, oth, nSpecTot):
        if(self.nSites != oth.nSites):
            return 0

        ##Make sorted species arrays
        s1, s2 = [0 for i in range(0, nSpecTot)], [0 for i in range(0, nSpecTot)]
        for i in range(0, self.nSites):
            s1[self.sites[i].spec] += 1
        for i in range(0, oth.nSites):
            s2[oth.sites[i].spec] += 1
        ##Check for identical species array
        for i in range(0, nSpecTot):
            if(s1[i] != s2[i]):
                return 0

        ##Make sorted distance arrays
        size = (self.nSites)*(self.nSites - 1) // 2
        d1, d2 = [None for i in range(0, size)], [None for i in range(0, size)]
        loc = 0
        for i in range(0, self.nSites - 1):
            for j in range(i + 1, self.nSites):
                d1[loc] = DistSqd_s(self.sites[i], self.sites[j])
                loc += 1
        d1.sort()
        loc = 0
        for i in range(0, oth.nSites - 1):
            for j in range(i + 1, oth.nSites):
                d2[loc] = DistSqd_s(oth.sites[i], oth.sites[j])
                loc += 1
        d2.sort()
        ##Check for identical distances
        for i in range(0, size):
            if(abs(d1[i] - d2[i]) > CART_EPS):
                return 0

        #Finially, if all of these checks pass, the two envs are identical
        return 1

    def AddSiteWithSpec(self, siteLine):
        self.nSites += 1
        s = site()
        s.InitWithSpec(lineString=siteLine)
        self.sites.append(s)

    def AddSiteWithoutSpec(self, siteLine):
        self.nSites += 1
        s = site()
        s.InitWithoutSpec(lineString=siteLine)
        self.sites.append(s)

    def SetEnergy(self, nrg):
        self.nrg = nrg

    def SetRep(self, rep):
        self.isRepped = rep

    #Write an env to a poscar file
    #Assumes that site species are initialized to strings representing atom types
    def ToPoscar(self, outLoc, A, comment="env printout"):
        ##Count elements
        elemCounts = dict()
        for sit in self.sites:
            elemCounts[sit.spec] = 0
        for sit in self.sites:
            elemCounts[sit.spec] += 1

        with open(outLoc, 'w') as outfile:
            ##Header
            outfile.write(str(comment) + "\n")
            for i in range(0, 3):
                for j in range(0, 3):
                    outfile.write(ffmt(A[i][j]) + ' ')
                outfile.write("\n")
            vas = "" ###keep this in order that keys are iterated: python makes no promises on the order
            ats = [] ###of a hash iteration as far as I know
            for ke, va in elemCounts.items():
                outfile.write(ke + ' ')
                ats.append(ke)
                vas += str(va) + ' '
            outfile.write("\n" + vas + "\nDi\n")
            ##Atom positions
            for at in ats:
                for sit in self.sites:
                    if(sit.spec == at):
                        for i in range(0, 3):
                            outfile.write(ffmt(sit.crdsD[i]) + ' ')
                        outfile.write("\n")
            outfile.close()

###
###  File Input
###
#Reads a lattice.tsam file and returns the real-space lattice vector matrix (A), a list of atom sites
#(not seperated into sublattices or anything crazy), and a mapping of species such that spec[1] = Ag
#means that species number 1 is mapped to Ag.
def ReadLatticeFile(infileLoc):
    A = [[None for j in range(0, 3)] for i in range(0, 3)]
    allSites = []
    map = []

    with open(infileLoc, 'r') as infile:
        cell, sublattice = 0, 0
        for lin in infile:
            if(cell + sublattice == 2):
                break

            #Lattice vector matrix
            if("BEGIN CELL" in lin):
                for i in range(0, 3):
                    line = infile.readline().split()
                    for j in range(0, 3):
                        A[i][j] = float(line[j])
                cell = 1
                continue

            #Species map, sites
            if("BEGIN SUBLATTICE" in lin):
                _ = infile.readline()
                map = [i for i in infile.readline().split()]

                nSublats = int(infile.readline())
                for i in range(0, nSublats):
                    _ = infile.readline()

                    nSites, stri = 0, [s for s in infile.readline().split()]
                    for n, s in enumerate(stri):
                        if(n%2 != 0):
                            nSites += int(s)
                    for j in range(0, nSites):
                        si = site()
                        si.InitWithoutSpec(infile.readline(), cCrds=False)
                        allSites.append(si)

                sublattice = 1
                continue

        infile.close()
    return A, allSites, map


#Reads the sparse-formatted env decomp file and returns a (dense!) double array of its entries
#i.e. in the sparse file, M[i][j] = 44, M[i][j+1] = 3 means that there are 3 instances of
#unique env number 44 in unique occupancy configuration i
#The returned matrix says that if M'[i][j] = 3, there are three instances of unique env number j in
#unique occupancy configuration i
def ReadDecompFile(infileLoc):
    ret = []
    with open(infileLoc, 'r') as infile:
        dim = [int(i) for i in infile.readline().split()]
        r, c = dim[0], dim[1]
        ret = [[0 for j in range(0, c)] for i in range(0, r)]

        for i in range(0, r):
            row = [int(k) for k in infile.readline().split()]
            for j in range(0, len(row), 2):
                ret[i][row[j]] = row[j + 1]

        infile.close()
    return ret

#Reads an environment printout file and returns a list of environments ordered
def ReadEnvFile(infileLoc):
    ret = []
    with open(infileLoc, 'r') as infile:
        nEnvs = int(infile.readline())
        for i in range(0, nEnvs):
            nSites = int(infile.readline())
            e = env()
            for j in range(0, nSites):
                e.AddSiteWithSpec(siteLine=infile.readline())
            _ = infile.readline() ##skip energy - it's NULL for now
            e.SetRep(rep=int(infile.readline()))
            ret.append(e)

        infile.close()
    return ret

#Reads a site-occupancy file and returns a list of site occupancies
def ReadOccFile(infileLoc):
    ret = []
    with open(infileLoc, 'r') as infile:
        dim = [int(i) for i in infile.readline().split()]
        r, c = dim[0], dim[1]
        ret = [[0 for j in range(0, c)] for i in range(0, r)]

        for i in range(0, r):
            ret[i] = [int(j) for j in infile.readline().split()]

        infile.close()
    return ret

###
###  File Output
###
#Writes a list of environments in a format that the simulated annealing program can deal with
#Assumes that we've already set the energies of the environments
def WriteEnvFile(outfileLoc, envs, prec=DEFAULT_PREC):
    with open(outfileLoc, 'w') as outfile:
        outfile.write("BEGIN ENVS\n")
        outfile.write(str(len(envs)) + "\n")
        for e in envs:
            outfile.write(str(e.nSites) + "\n")
            for s in e.sites:
                outfile.write(str(s.spec))
                for i in range(0, 3):
                    outfile.write(' ' + ffmt(s.crdsD[i], prec))
                for i in range(0, 3):
                    outfile.write(' ' + ffmt(s.crdsC[i], prec))
                outfile.write("\n")
            outfile.write(ffmt(e.nrg, prec) + "\n")
            outfile.write(str(e.isRepped) + "\n")
        outfile.write("END ENVS\n")
        outfile.close()


#Makes a poscar from a list of sites, cell vectors, occupancy vector and a mapping of ints -> species
#Option to exclude ceritan elements (Enter as a LIST i.e. ["Va", "Vac"] to exclude vacancies regardless of
#whether they're represented as 'Va' (a bad idea, BTW), or 'Vac').
def WritePoscar_q(outfileLoc, A, sites, occs, map, exclude=[], comment="default comment"):
    ##Setup: make sure stuff is in the right order
    orderedInts, orderedSpecs, orderedSpecCnts = [], [], []
    for n, m in enumerate(map):
        if m not in exclude:
            orderedInts.append(n)
            orderedSpecs.append(m)
            orderedSpecCnts.append(0)

    orderedSites = []
    for i in range(0, len(orderedInts)):
        for j in range(0, len(occs)):
            if(orderedInts[i] == occs[j]):
                orderedSites.append(sites[j])
                orderedSpecCnts[i] += 1

    ##Make the acutal POSCAR
    with open(outfileLoc, 'w') as outfile:
        ###Header
        outfile.write(comment + "\n1.0\n")
        for i in range(0, 3):
            for j in range(0, 3):
                outfile.write(ffmt(A[i][j]) + ' ')
            outfile.write("\n")
        ###Atoms
        for spec in orderedSpecs:
            outfile.write(spec + ' ')
        outfile.write("\n")
        for specC in orderedSpecCnts:
            outfile.write(str(specC) + ' ')
        outfile.write("\nDi\n")
        ###Atom positions
        for sit in orderedSites:
            for i in range(0, 3):
                outfile.write(ffmt(sit.crdsD[i]) + ' ')
            outfile.write("\n")
        outfile.close()
