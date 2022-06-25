#if 0
//Basic Simulated Annealing Algo
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>

#include "Decs.h"
#include "Constants.h"

//min only defined in windows VS somehow
#define min(X,Y) (((X) < (Y)) ? (X) : (Y))

char INFILE_LATT[LINESIZE] = "lattice.tsam\0";
char INFILE_PRMS[LINESIZE] = "runparams.tsam\0";
char INFILE_ENVS[LINESIZE] = "envs.tsam\0";

int main(int argc, char *argv[]){
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
		//Infile (environments) name
		if((argv[i][0] == '-') && (argv[i][1] == 'e' || argv[i][1] == 'E')){
			for(ushort j = 0; j < LINESIZE; ++j){
				if(argv[i + 1][j] == ' '){
					break;
				}
				INFILE_ENVS[j] = argv[i + 1][j];
			}
		}

		//Overwritten cutoff radius
		if((argv[i][0] == '-') && (argv[i][1] == 'r' || argv[i][1] == 'R')){
			ovwRCut = atof(argv[i + 1]);
		}

		//Overwritten annealing parameters
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

		//Overwritten swapping probabilities (must be normalized to sum to 1)
		if((argv[i][0] == '-') && (argv[i][1] == 's' || argv[i][1] == 'S')){
			int nEntries = atol(argv[i + 1]);
			ovwSwapProbs = malloc(nEntries * sizeof(double));
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
	if(ovwNStep > 0u){
		nSteps = ovwNStep; startStep = ovwSStep;
		startNrg = ovwSNrg; endNrg = ovwENrg;
		startCool = ovwSCool; endCool = ovwECool;
		alpha = ovwAlpha;
		rSeed = ovwRSeed;
	}
	printf("ok\n");
	
	printf("LOG: reading environments from file %s ... ", INFILE_ENVS);
	uint nFitEnvs; env** fitEnvs;
	ushort* fixedSpecArr = malloc(nSpecies*sizeof(ushort));
	for(ushort i = 0; i < nSpecies; ++i) fixedSpecArr[i] = i;
	ReadEnvFile(INFILE_ENVS, &nFitEnvs, &fitEnvs,
				nSpecies, fixedSpecArr);
	printf("ok\n");

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
	for(ushort i = 0; i < nSublats; ++i) absPs[i] = sublats[i]->swapProb;
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

	//Holds info on any unknown structures (needs to be before any gotos)
	ushort nUnkEnvs; ushort* unkEnvs = calloc(nSites, sizeof(ushort));
	//Ditto for the best site printout
	ushort* bestSites = malloc(nSites * sizeof(ushort));

	//Predict the energy of init state
	double e0 = CalcNrg(nFitEnvs, fitEnvs, nSites, &envBase,
						&nUnkEnvs, &unkEnvs, nSpecies);
	if(nUnkEnvs){
		UnknownEnvHandler(nSites, allSites, 0, nUnkEnvs, unkEnvs, 
						  nSites, envBase); printf("\n");
		goto Clean_Exit;
	}
	e0 /= nSites;
	printf("JOB: energy of initial configuration e0 = %f eV/site\n", e0);

	//Simulated annealing------------------------------------------------------
	printf("OPT: n%7c   e(n)%7c   e(n-1)%5c   e(opt)%5c   tau%8c   p(swap)\n", 
		   ' ', ' ', ' ', ' ', ' ');
	printf("OPT: ------------------------------------------------------------"
		   "--------------\n");
	
	ushort swapSublat = 0; ushort swapSpecInd1 = 0; ushort swapSpecInd2 = 0;
	double eOpt = e0;
	double eDiff;
	double tau;
	double prob;
	for(uint i = startStep; i < nSteps; ++i){
		///Trial step
		UniformRandSwaps(nSublats, &sublats, mapLo, mapHi, N_NBR_SWAPS, 
						 &swapSublat, &swapSpecInd1, &swapSpecInd2);
		double eT = CalcNrg(nFitEnvs, fitEnvs, nSites, &envBase, 
							&nUnkEnvs, &unkEnvs, nSpecies);
		if(nUnkEnvs){
			UnknownEnvHandler(nSites, allSites, i, nUnkEnvs, unkEnvs,
							  nSites, envBase); 
			printf("\n");
			goto Clean_Exit;
		}
		eT /= nSites;

		///Decide if we should accept new structure
		eDiff = eT - e0;
		tau = CalcTau(i, nSteps, startNrg, endNrg, startCool, endCool, alpha);
		prob = min(exp(-eDiff/tau), 1.0);
		printf("OPT: %08u   %011.7f   %011.7f   %011.7f   %011.6f   %06.5f\n",
			   i + 1u, eT, e0, eOpt, tau, prob);
		///Always swap for lower energy state,
		///sometimes swap for higher energy state
		if((((double)rand()/(double)RAND_MAX) < prob) || (eDiff < 0.0)){
			e0 = eT;
		}
		else{
			///If we're here, don't swap - move back to orig struct
			ushort tmp = *(sublats[swapSublat]->sites[swapSpecInd1]->species);
			*(sublats[swapSublat]->sites[swapSpecInd1]->species) = 
						*(sublats[swapSublat]->sites[swapSpecInd2]->species);
			*(sublats[swapSublat]->sites[swapSpecInd2]->species) = tmp;
		}
		
		///Update optimals
		if(e0 < eOpt){ 
			eOpt = eT;
			for(ushort j = 0; j < nSites; ++j){
				bestSites[j] = *allSites[j]->species;
			}
		}
	}

	//Print out the optimal energy and configuration found
	printf("JOB: energy of final configuration eOpt = %f eV/site\n", eOpt);
	printf("JOB: BEGIN OPT SITE PRINTOUT\n");
	printf("%u %u\n", 1u, nSites);
	for(ushort i = 0; i < nSites; ++i) printf("%u ", bestSites[i]);
	printf("\nJOB: END OPT SITE PRINTOUT\n");

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
	for(uint i = 0; i < nFitEnvs; ++i){
		for(ushort j = 0; j < fitEnvs[i]->nSites; ++j){
			free(fitEnvs[i]->sites[j]);
		}
		free(fitEnvs[i]->sortedSpecies);
		free(fitEnvs[i]->sortedDists);
		free(fitEnvs[i]);
	}
	free(fitEnvs);
	for(ushort i = 0; i < nSublats; ++i){
		free(sublats[i]->sites);
		free(sublats[i]);
	}
	free(sublats);
	free(fixedSpecArr);
	free(initState);
	free(mapLo);
	free(mapHi);
	free(unkEnvs);
	free(bestSites);

	return 0;
}
#endif