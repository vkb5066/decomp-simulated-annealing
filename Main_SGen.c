#if 1
//Uniform Random Structure Generation
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "Decs.h"
#include "Constants.h"


char INFILE_LATT[LINESIZE] = "lattice.tsam\0";
char INFILE_PRMS[LINESIZE] = "runparams.tsam\0";
char INFILE_WARM[LINESIZE] = "warmstart.tsam\0";


int main(int argc, char *argv[]){
	ushort warmStart = 0;
	double ovwRCut = -1.0;
	uint ovwNStep = 0; uint ovwSStep; double ovwSNrg; double ovwENrg;
	double ovwSCool; double ovwECool; double ovwAlpha; int ovwRSeed;
	double* ovwSwapProbs = NULL;
	//CLIs---------------------------------------------------------------------
	for (ushort i = 0; i < argc; ++i){
		//Infile (lattice) name
		if((argv[i][0] == '-') && (argv[i][1] == 'l' || argv[i][1] == 'L')){
			for(ushort j = 0; j < LINESIZE; ++j){
				if(argv[i + 1][j] == ' '){
					break;
				}
				INFILE_LATT[j] = argv[i + 1][j];
			}
		}
		//Infile (parameters) name
		if((argv[i][0] == '-') && (argv[i][1] == 'p' || argv[i][1] == 'P')){
			for(ushort j = 0; j < LINESIZE; ++j){
				if(argv[i + 1][j] == ' '){
					break;
				}
				INFILE_PRMS[j] = argv[i + 1][j];
			}
		}
		//Infile (warm restart) name, flag: envs
		if((argv[i][0] == '-') && (argv[i][1] == 'w' || argv[i][1] == 'W')){
			warmStart = (ushort)1;
			for(ushort j = 0; j < LINESIZE; ++j){
				if(argv[i + 1][j] == ' '){
					break;
				}
				INFILE_WARM[j] = argv[i + 1][j];
			}
		}

		//Overwritten cutoff radius
		if((argv[i][0] == '-') && (argv[i][1] == 'r' || argv[i][1] == 'R')){
			ovwRCut = atof(argv[i + 1]);
		}

		//Overwritten annealing parameters (most of these are unused - keep
		//for compatability w/ simulated annealing routine)
		if((argv[i][0] == '-') && (argv[i][1] == 'a' || argv[i][1] == 'A')){
			ovwNStep = (uint)atol(argv[i + 1]);
			ovwSStep = (uint)atol(argv[i + 2]);
			ovwSNrg = atof(argv[i + 3]);
			ovwENrg = atof(argv[i + 4]);
			ovwSCool = atof(argv[i + 5]);
			ovwECool = atof(argv[i + 6]);
			ovwAlpha = atof(argv[i + 7]);
			ovwRSeed = atol(argv[i + 8]);
		}

		//Overwritten number of gen steps (easier than the -a command)
		if((argv[i][0] == '-') && (argv[i][1] == 'n' || argv[i][1] == 'N')){
			ovwNStep = (uint)atol(argv[i + 1]);
		}

		//Overwritten swapping probabilities (must be normalized to sum to 1)
		if((argv[i][0] == '-') && (argv[i][1] == 's' || argv[i][1] == 'S')){
			int nEntries = atol(argv[i + 1]);
			ovwSwapProbs = malloc(nEntries*sizeof(double));
			for(ushort j = 0; j < nEntries; ++j){
				ovwSwapProbs[j] = atof(argv[i + 2 + j]);
			}
		}
	}
	//-------------------------------------------------------------------------


	//Setup--------------------------------------------------------------------
	//Read input files
	printf("LOG: reading lattice from file %s ... ", INFILE_LATT);
	double A[3][3];
	ushort nSpecies;
	ushort nSites; site** allSites;
	ushort nSublats; sublat** sublats; 
	ReadLatticeFile(INFILE_LATT, &A, &nSpecies, &nSites, &allSites, 
					&nSublats, &sublats);
	printf("ok\n");

	printf("LOG: reading run params from file %s ... ", INFILE_PRMS); 
	double rCut;
	uint nSteps; uint startStep; double startNrg; double endNrg;
	double startCool; double endCool; double alpha; int rSeed;
	ushort* initState = malloc(nSites*sizeof(ushort));
	ushort readInitConfig;
	readInitConfig = ReadParamFile(INFILE_PRMS, &rCut, &nSteps, &startStep,
		&startNrg, &endNrg, &startCool, &endCool,
		&alpha, &rSeed, &initState);
	if(ovwRCut > 0.0) rCut = ovwRCut;
	if(ovwNStep > 0u) nSteps = ovwNStep;
	printf("ok\n");
	
	ushort* fixedSpecArr = malloc(nSpecies*sizeof(ushort));
	for (ushort i = 0; i < nSpecies; ++i) fixedSpecArr[i] = i;

	//Fetch chemical env geoms (and set cartesian coordinates for imgs)
	printf("LOG: fetching chemical environment geometries ... ");
	for(ushort i = 0; i < nSites; ++i) MoveToCell(&allSites[i]);
	ushort imgArr[3]; GiveNImgs(&imgArr, rCut, A);
	env** envBase;
	SetSiteGeoms(nSites, &envBase, nSites, allSites, imgArr, rCut, A, 
				 nSpecies);
	printf("ok\n");


	//Initialization-----------------------------------------------------------
	if(rSeed == -1) rSeed = time(NULL);
	srand((uint)rSeed);
	printf("JOB: rng seed = %u\n", (uint)rSeed);

	//Mapping of probabilities for uniform random choices of sublattices
	double* mapLo; double* mapHi;
	double* absPs = malloc(nSublats * sizeof(double));
	for (ushort i = 0; i < nSublats; ++i) absPs[i] = sublats[i]->swapProb;
	UniformMap(&mapLo, &mapHi, nSublats, absPs);
	free(absPs);

	//Set initial configuration state if provided.  Otherwise randomize it
	printf("JOB: setting initial configuration ");
	if(readInitConfig) printf("from file %s ... ", INFILE_PRMS);
	else printf("randomly ... ");
	if(readInitConfig) SetConfig(nSites, &allSites, initState);
	else NonUniformRandSwaps(nSublats, &sublats);
	printf("ok\n");

	//Show absolute probabilities of sublattice species swaps
	for(ushort i = 0; i < nSublats; ++i){
		if(ovwSwapProbs) sublats[i]->swapProb = ovwSwapProbs[i];
		printf("JOB: p(s%u) = %f\n", i, sublats[i]->swapProb);
	}
	free(ovwSwapProbs);
	//Avoid possible cataclysmic disaster
	if(N_REALLOC < 2*nSites + 1){
		printf("ERR: N_REALLOC < 2(number of sites) + 1\n");
		exit(1);
	}

	//Structure generation-----------------------------------------------------
	//If a random site config vector results in a unique chem env vector, keep 
	//both of them

	///Warm start: re-initialize stuff if necessary
	////unique environment array
	uint nUnqEnvC; uint nUnqEnvM; env** unqEnvs;
	uint nUnreppedEnvs;
	////unique chem env decompositions
	uint nUnqDecC; uint nUnqDecM; uint** unqDecs;
	///stuff for unique site "decompositions" - in line with unq env decomps
	ushort** unqSites;
	
	if(warmStart){
		printf("LOG: warm start - reading from file %s ... ",
			   INFILE_WARM);
		ReadWarmStartFile(INFILE_WARM, &nUnqEnvC, &nUnqEnvM, &unqEnvs,
						  &nUnqDecC, &nUnqDecM, &unqDecs, &unqSites,
						  &nUnreppedEnvs, nSpecies, fixedSpecArr, nSites);
		printf("ok\n");
	}
	else{	
		printf("LOG: cold start\n");

		nUnqEnvC = 0u; ////(c)urrent # unique envs
		nUnqEnvM = N_REALLOC; ////(m)ax # unique envs that I have space for
		unqEnvs = malloc(N_REALLOC*sizeof(env*));
		nUnreppedEnvs = 0u;

		nUnqDecC = 0u; ////(c)urrent # unique decomps
		nUnqDecM = N_REALLOC; ////(m)ax # unique decomps that I have space for
		unqDecs = malloc(N_REALLOC*sizeof(uint*));
	
		unqSites = malloc(N_REALLOC*sizeof(ushort*));
	}
	ushort* currDecArr = malloc(nUnqEnvM*sizeof(ushort));

	///Continue to generate and eval random site permutations
	printf("GEN: n(tot)   n(unqS)  n(unqE)  n(lunqE) n(unrE)\n");
	ushort swapSublat = 0; ushort swapSpecInd1 = 0; ushort swapSpecInd2 = 0;
	uint slue = 0; ///since last unique env
	for(uint n = 0; n < nSteps; ++n){
		////set the current cell's decomposition...
		memset(currDecArr, (ushort)0, nUnqEnvM*sizeof(ushort));
		for(ushort i = 0; i < nSites; ++i){
			ushort toAdd = 1;
			SetSpecArr_e(&envBase[i], nSpecies);
			SetDistArr_e(&envBase[i]);

			for(uint j = 0; j < nUnqEnvC; ++j){
				if(Cmpr_e(envBase[i], unqEnvs[j], nSpecies)){
					toAdd = (ushort)0;
					currDecArr[j]++;
					break;
				}
			}

			/////add a copy of this environment if it is brand new
			if(toAdd){
				slue = 0;
				currDecArr[nUnqEnvC]++;

				env* newEnv = malloc(sizeof(env));
				newEnv->nSites = envBase[i]->nSites;
				newEnv->sites = malloc((newEnv->nSites)*sizeof(site));
				for(ushort j = 0; j < newEnv->nSites; ++j){
					site* newSite = malloc(sizeof(site));
					newSite->species = fixedSpecArr + 
									   *(envBase[i]->sites[j]->species);
					for(ushort m = 0; m < 3; ++m){
						newSite->crdsD[m] = envBase[i]->sites[j]->crdsD[m];
						newSite->crdsC[m] = envBase[i]->sites[j]->crdsC[m];
					}
					newSite->self = envBase[i]->sites[j]->self;
					newEnv->sites[j] = newSite;
				}
				newEnv->sortedSpecies = malloc(nSpecies*sizeof(ushort));
				SetSpecArr_e(&newEnv, nSpecies);
				newEnv->nDists = envBase[i]->nDists;
				newEnv->sortedDists = malloc((newEnv->nDists)*sizeof(double));
				SetDistArr_e(&newEnv);

				newEnv->repped = (ushort)0;
				nUnreppedEnvs++;

				unqEnvs[nUnqEnvC] = newEnv;
				nUnqEnvC++;
			}
		}

		if(nUnreppedEnvs){
			ushort isUniqueDecomp = 1;
			for(uint i = 0; i < nUnqDecC; ++i){
				if(DecompCmpr_fs(nSites, currDecArr, unqDecs[i])){
					isUniqueDecomp = (ushort)0;
					break;
				}
			}

			////store sparse version of decomp + update repped envs if unique
			////and important
			if(isUniqueDecomp){
				unqDecs[nUnqDecC] = malloc(2*nSites*sizeof(uint));
				DecompConv(nSites, currDecArr, &(unqDecs[nUnqDecC]));
				ushort accum = 0;
				for(ushort i = 0; accum < nSites; i += (ushort)2){
					if(!(unqEnvs[unqDecs[nUnqDecC][i]]->repped)){
						unqEnvs[unqDecs[nUnqDecC][i]]->repped = (ushort)1;
						unqSites[nUnqDecC] = malloc(nSites*sizeof(ushort));
						for(ushort i = 0; i < nSites; ++i){
							unqSites[nUnqDecC][i] = *allSites[i]->species;
						}
						nUnreppedEnvs--;
						nUnqDecC++;

						goto Keep; //////only want one rep per unique decomp
					}
					accum += unqDecs[nUnqDecC][i + (ushort)1];
				}
				free(unqDecs[nUnqDecC]);
				Keep: NOP
			}
		}

		////allocate more space if necessary
		if(nUnqEnvM - nUnqEnvC <= nSites){ /////wont overflow: max >= current
			nUnqEnvM += N_REALLOC;
			unqEnvs = realloc(unqEnvs, nUnqEnvM*sizeof(env*));
			currDecArr = realloc(currDecArr, nUnqEnvM*sizeof(ushort));
		}
		if(nUnqDecC == nUnqDecM){
			nUnqDecM += N_REALLOC;
			unqDecs = realloc(unqDecs, nUnqDecM*sizeof(uint*));
			unqSites = realloc(unqSites, nUnqDecM*sizeof(ushort*));
		}
		
		////make new random permutation, print out info
		UniformRandSwaps(nSublats, &sublats, mapLo, mapHi, N_NBR_SWAPS,
						 &swapSublat, &swapSpecInd1, &swapSpecInd2);
		printf("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b"
			   "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b");
		printf("GEN: %08u %08u %08u %08u %08u", n, nUnqDecC,
			   nUnqEnvC, slue, nUnreppedEnvs);
		slue++;
	}
	printf("\n");


	//At the end, print out all:
	//Chem envs
	printf("JOB: BEGIN UNIQUE ENV PRINTOUT\n");
	EnvPrintout(nUnqEnvC, unqEnvs, (ushort)0, (ushort)1);
	printf("JOB: END UNIQUE ENV PRINTOUT\n");
	//Env decompositions
	printf("JOB: BEGIN UNIQUE ENV DECOMP PRINTOUT\n");
	DecPrintout(nUnqDecC, nUnqEnvC, nSites, unqDecs);
	printf("JOB: END UNIQUE ENV DECOMP PRINTOUT\n");
	//Site decompositions
	printf("JOB: BEGIN UNIQUE SITE PRINTOUT\n");
	SitePrintout(nUnqDecC, nSites, unqSites);
	printf("JOB: END UNIQUE SITE PRINTOUT\n");


	//Clean up-----------------------------------------------------------------
	Clean_Exit:
	printf("LOG: ending program\n");
	for(ushort i = 0; i < nSites; ++i){
		free(allSites[i]->species);
		free(allSites[i]);
		for(ushort j = 0; j < envBase[i]->nSites; ++j){
			free(envBase[i]->sites[j]);
		}
		free(envBase[i]->sortedSpecies);
		free(envBase[i]->sortedDists);
		free(envBase[i]);
	}
	free(allSites); 
	free(envBase);
	free(currDecArr);
	for(uint i = 0; i < nUnqEnvC; ++i){
		for(ushort j = 0; j < unqEnvs[i]->nSites; ++j){
			free(unqEnvs[i]->sites[j]);
		}
		free(unqEnvs[i]->sortedSpecies);
		free(unqEnvs[i]->sortedDists);
		free(unqEnvs[i]);
	}
	free(unqEnvs);
	for(uint i = 0; i < nUnqDecC; ++i){
		free(unqDecs[i]);
		free(unqSites[i]);
	}
	free(unqDecs);
	free(unqSites);
	for(ushort i = 0; i < nSublats; ++i){
		free(sublats[i]->sites);
		free(sublats[i]);
	}
	free(sublats);
	free(fixedSpecArr);
	free(initState);
	free(mapLo);
	free(mapHi);

	return 0;
}
#endif