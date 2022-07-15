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
	uvlong ovwNStep = 0llu; uvlong ovwSStep = 0llu; 
	double ovwSNrg = -1.0; int ovwRSeed;
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
			ovwNStep = strtoull(argv[i + 1], NULL, BASE);
			ovwSStep = strtoull(argv[i + 2], NULL, BASE);
			ovwSNrg = atof(argv[i + 3]);
			ovwRSeed = atol(argv[i + 4]);
		}

		//Overwritten number of gen steps (easier than the -a command)
		if((argv[i][0] == '-') && (argv[i][1] == 'n' || argv[i][1] == 'N')){
			ovwNStep = strtoull(argv[i + 1], NULL, BASE);
		}

		//Overwritten starting thermal energy (easier than the -a command)
		if((argv[i][0] == '-') && (argv[i][1] == 't' || argv[i][1] == 'T')){
			ovwSNrg = atof(argv[i + 1]);
		}

		//Overwritten swapping probabilities (must be normalized to sum to 1)
		if((argv[i][0] == '-') && (argv[i][1] == 's' || argv[i][1] == 'S')){
			ushort nEntries = (ushort)atol(argv[i + 1]);
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
	ushort nSites; site** allSites; site** allSitesTrial;
	ushort nSublats; sublat** sublats; sublat** sublatsTrial;
	ReadLatticeFile(INFILE_LATT, &A, &nSpecies, &nSites, &allSites, 
					&nSublats, &sublats);
	ReadLatticeFile(INFILE_LATT, &A, &nSpecies, &nSites, &allSitesTrial,
					&nSublats, &sublatsTrial);
	printf("ok\n");

	printf("LOG: reading run params from file %s ... ", INFILE_PRMS);
	double rCut;
	uvlong nSteps; uvlong startStep; double startNrg; int rSeed;
	ushort* initState = malloc(nSites*sizeof(ushort));
	ushort readInitConfig;
	readInitConfig = ReadParamFile(INFILE_PRMS, &rCut, &nSteps, &startStep,
								   &startNrg, &rSeed, &initState);
	if(ovwRCut > 0.0) rCut = ovwRCut;
	if(ovwNStep > 0u) nSteps = ovwNStep;
	if(ovwSNrg > 0.0) startNrg = ovwSNrg;
	if(ovwSStep != 0u){
		nSteps = ovwNStep; startStep = ovwSStep;
		startNrg = ovwSNrg;
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
	env** envBase; env** envBaseTrial;
	SetSiteGeoms(nSites, &envBase, nSites, &allSites, imgArr, rCut, A, 
				 nSpecies);
	SetSiteGeoms(nSites, &envBaseTrial, nSites, &allSitesTrial, imgArr, rCut, 
				 A, nSpecies);
	printf("ok\n");

	//Make the hash-tree structure
	printf("LOG: creating hash-tree ... ");
	hbkt table[HASH_TABLE_LEN]; InitTable_h(&table);
	ushort colls = 0; uint successes = 0; env* tmp;
	for(uint i = 0; i < nFitEnvs; ++i)
		tmp = AddEnv(&table, nSpecies, fitEnvs[i], &colls, &successes);
	///make sure my algorithm doesn't suck
	if(successes != nFitEnvs){
		printf("\nERR: unable to add all fit envs to hash-tree\n");
		return 1;
	}
	printf("ok (%u collisions)\n", colls);


	//Initialization-----------------------------------------------------------
	if(rSeed == -1) rSeed = time(NULL);
	srand((uint)rSeed);
	printf("JOB: rng seed = %u\n", (uint)rSeed);

	//Set initial configuration state if provided.  Otherwise randomize it
	printf("JOB: setting initial configuration ");
	if(readInitConfig) printf("from file %s ... ", INFILE_PRMS);
	else printf("randomly ... ");
	if(readInitConfig) SetConfig(nSites, &allSites, initState);
	else NonUniformRandSwaps(nSublats, &sublats);
	DeepcopyState(nSites, allSites, &allSitesTrial);
	printf("ok\n");

	//Show absolute probabilities of sublattice species swaps
	//also overwrite default swapping probabilities if necessary
	for(ushort i = 0; i < nSublats; ++i){
		if(ovwSwapProbs){
			sublats[i]->swapProb = ovwSwapProbs[i];
			sublatsTrial[i]->swapProb = ovwSwapProbs[i];
		}
		printf("JOB: p(s%u) = %f\n", i, sublats[i]->swapProb);
	}
	free(ovwSwapProbs);

	//Mapping of probabilities for uniform random choices of sublattices
	double* mapLo; double* mapHi;
	double* absPs = malloc(nSublats*sizeof(double));
	for(ushort i = 0; i < nSublats; ++i) absPs[i] = sublats[i]->swapProb;
	UniformMap(&mapLo, &mapHi, nSublats, absPs);
	free(absPs);

	//Holds info on any unknown structures (needs to be before any gotos)
	ushort nUnkEnvs; ushort* unkEnvs = calloc(nSites, sizeof(ushort));
	//Ditto for the best site printout
	ushort* bestSites = malloc(nSites*sizeof(ushort));

	//Mark ALL initial envs as needing to be calculated + set their
	//energies to zero
	for(ushort i = 0; i < nSites; ++i){ 
		envBase[i]->recalc = (ushort)1;
		envBase[i]->nrg = 0.0;
		envBaseTrial[i]->recalc = (ushort)1;
		envBaseTrial[i]->nrg = 0.0;
	}

	//Array of indices that mark which envs from the base envs have been 
	//recalculated - needed to avoid looping over all sites after deciding
	//whether to accept or reject a move
	ushort nRecalcInds = 0;
	ushort *recalcInds = malloc(nSites*sizeof(ushort));

	//Predict the energy of init state
	double eCurr = 0.0;
	UpdateNrg(&eCurr, nSites, table, nSites, &envBase,
			  &nUnkEnvs, &unkEnvs, nSpecies, &nRecalcInds, &recalcInds);
	eCurr = 0.0;
	UpdateNrg(&eCurr, nSites, table, nSites, &envBaseTrial,
			  &nUnkEnvs, &unkEnvs, nSpecies, &nRecalcInds, &recalcInds);
#if !(ASSIGN_UNK_ENV_NRG)
	if(nUnkEnvs){
		UnknownEnvHandler(nSites, allSites, 0, nUnkEnvs, unkEnvs,
						  nSites, envBase); printf("\n");
		goto Clean_Exit;
	}
#endif
	printf("JOB: energy of initial configuration e0 = %f eV/site\n", eCurr);
	//and initialize the best sites to the initial sites
	for(ushort i = 0; i < nSites; ++i){
		bestSites[i] = *allSites[i]->species;
	}	
	

	//Simulated annealing------------------------------------------------------
	printf("OPT: n%10c  e(n)%7c  e(n-1)%5c  e(opt)%5c  tau%8c  noise\n", 
		   ' ', ' ', ' ', ' ', ' ');
	printf("OPT: ------------------------------------------------------------"
		   "------------\n");
	
	ushort swpS = 0; ushort swpI1 = 0; ushort swpI2 = 0;
	double eOpt = eCurr; double eTrial = eCurr;
	double eDiff;
	double tau;
	double prob;
	double noise = 1.0;
	ushort t;
	for(uvlong i = startStep; i < nSteps; ++i){

		///Trial step
		UniformRandSwaps(nSublats, &sublatsTrial, mapLo, mapHi, 1u, 
						 &swpS, &swpI1, &swpI2);
		MarkEnvsToRecalc(sublatsTrial, swpS, swpI1, swpI2);
		UpdateNrg(&eTrial, nSites, table, nSites, &envBaseTrial, 
				  &nUnkEnvs, &unkEnvs, nSpecies, &nRecalcInds, &recalcInds);
#if !(ASSIGN_UNK_ENV_NRG)
		if(nUnkEnvs){
			UnknownEnvHandler(nSites, allSitesTrial, i, nUnkEnvs, unkEnvs,
							  nSites, envBaseTrial); printf("\n");
			goto Clean_Exit;
		}
#endif
		
		///Decide if we should accept new structure
		tau = CalcTau(i, nSteps, startNrg);
		tau *= noise;
		eDiff = eTrial - eCurr;
		prob = min(exp(-eDiff/tau), 1.0);

		//Print, if necessary
#if ANN_PRINTSTYLE
		if(i%ANN_PRINTEVERY == 0llu){
#if (ANN_PRINTSTYLE == 1)
			printf("OPT: %011llu  %011.7f  %011.7f  %011.7f  %011.9f  %07.5f"
				   "\n", i + 1llu, eTrial, eCurr, eOpt, tau, noise);
#endif
#if (ANN_PRINTSTYLE == 2)
			printf("OPT: %07.3f%% completed\r", (float)i/(float)nSteps*100.0f);
			fflush(stdout);
#endif
		}
#endif

		///Always swap for lower energy state,
		///sometimes swap for higher energy state
		if((((double)rand()/(double)RAND_MAX) < prob) || (eDiff < 0.0)){
			////Case: accept move
			eCurr = eTrial;
			*(sublats[swpS]->sites[swpI1]->species) = 
								*(sublatsTrial[swpS]->sites[swpI1]->species);
			*(sublats[swpS]->sites[swpI2]->species) = 
								*(sublatsTrial[swpS]->sites[swpI2]->species);
			////reset the envs
			for(ushort j = 0; j < nRecalcInds; ++j){
				t = recalcInds[j];
				envBase[t]->recalc = (ushort)0;
				SwapCmpArrs_e(envBaseTrial[t], &envBase[t], nSpecies);
				envBase[t]->nrg = envBaseTrial[t]->nrg;
			}	
		}
		else{
			////Case: reject move
			eTrial = eCurr;
			*(sublatsTrial[swpS]->sites[swpI1]->species) = 
									*(sublats[swpS]->sites[swpI1]->species);
			*(sublatsTrial[swpS]->sites[swpI2]->species) = 
									*(sublats[swpS]->sites[swpI2]->species);
			////reset the envs
			for(ushort j = 0; j < nRecalcInds; ++j){
				t = recalcInds[j];
				envBaseTrial[t]->recalc = (ushort)0;
				SwapCmpArrs_e(envBase[t], &envBaseTrial[t], nSpecies);
				envBaseTrial[t]->nrg = envBase[t]->nrg;
			}
		}
		
		///Update optimals, compute next iteration noise
		if(eCurr < eOpt){ 
			eOpt = eCurr;
			for(ushort j = 0; j < nSites; ++j){
				bestSites[j] = *allSites[j]->species;
			}
		}
		noise = 1.0 - (eCurr - eOpt)/eCurr; ///minus bc e0 is negative
	}
	printf("OPT: %011llu  %011.7f  %011.7f  %011.7f  %011.9f  %07.5f"
		   "\n", nSteps, eTrial, eCurr, eOpt, tau, noise);
	

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
		///current
		free(allSites[i]->species);
		free(allSites[i]->parentEnvs); ///free the actual envs later
		free(allSites[i]);
		for(ushort j = 0; j < envBase[i]->nSites; ++j){
			free(envBase[i]->sites[j]);
		}
		free(envBase[i]->sites);
		free(envBase[i]->sortedSpecies);
		free(envBase[i]->sortedDists);
		free(envBase[i]);
		///trial
		free(allSitesTrial[i]->species);
		free(allSitesTrial[i]->parentEnvs); ///free the actual envs later
		free(allSitesTrial[i]);
		for(ushort j = 0; j < envBaseTrial[i]->nSites; ++j){
			free(envBaseTrial[i]->sites[j]);
		}
		free(envBaseTrial[i]->sites);
		free(envBaseTrial[i]->sortedSpecies);
		free(envBaseTrial[i]->sortedDists);
		free(envBaseTrial[i]);
	}
	free(allSites); 
	free(envBase);
	free(allSitesTrial);
	free(envBaseTrial);
	for(uint i = 0; i < nFitEnvs; ++i){
		for(ushort j = 0; j < fitEnvs[i]->nSites; ++j){
			free(fitEnvs[i]->sites[j]);
		}
		free(fitEnvs[i]->sites);
		free(fitEnvs[i]->sortedSpecies);
		free(fitEnvs[i]->sortedDists);
		free(fitEnvs[i]);
	}
	free(fitEnvs);
	for(ushort i = 0; i < nSublats; ++i){
		///current
		free(sublats[i]->sites);
		free(sublats[i]);
		///trial
		free(sublatsTrial[i]->sites);
		free(sublatsTrial[i]);
	}
	free(sublats);
	free(sublatsTrial);
	free(fixedSpecArr);
	free(initState);
	DeallocTable_h(&table);
	free(mapLo);
	free(mapHi);
	free(unkEnvs);
	free(bestSites);
	free(recalcInds);

	return 0;
}
#endif