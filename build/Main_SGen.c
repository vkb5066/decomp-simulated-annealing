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
	uvlong ovwNStep = 0llu; uvlong ovwSStep;
	double ovwSNrg = -1.0; int ovwRSeed;
	ushort ovwNDecsPerEnv = 1;
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
			ovwNStep = strtoull(argv[i + 1], NULL, BASE);
			ovwSStep = strtoull(argv[i + 2], NULL, BASE);
			ovwSNrg = atof(argv[i + 3]);
			ovwRSeed = atol(argv[i + 4]);
		}

		//Overwritten number of gen steps (easier than the -a command)
		if((argv[i][0] == '-') && (argv[i][1] == 'n' || argv[i][1] == 'N')){
			ovwNStep = strtoul(argv[i + 1], NULL, BASE);
		}

		//Overwritten number of env decs per unique env 
		if((argv[i][0] == '-') && (argv[i][1] == 'd' || argv[i][1] == 'D')){
			ovwNDecsPerEnv = (ushort)atol(argv[i + 1]);
		}

		//Overwritten starting thermal energy (easier than the -a command)
		if((argv[i][0] == '-') && (argv[i][1] == 't' || argv[i][1] == 'T')){
			ovwSNrg = atof(argv[i + 1]);
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
	uvlong nSteps; uvlong startStep; double startNrg; int rSeed;
	ushort* initState = malloc(nSites*sizeof(ushort));
	ushort readInitConfig;
	readInitConfig = ReadParamFile(INFILE_PRMS, &rCut, &nSteps, &startStep,
								   &startNrg, &rSeed, &initState);
	ushort nDecsPerEnv = (ushort)startStep;
	if(ovwRCut > 0.0) rCut = ovwRCut;
	if(ovwNStep > 0llu) nSteps = ovwNStep;
	if(ovwNDecsPerEnv > (ushort)1) nDecsPerEnv = ovwNDecsPerEnv;
	if(ovwSNrg > 0.0) startNrg = ovwSNrg;
	printf("ok\n");
	
	ushort* fixedSpecArr = malloc(nSpecies*sizeof(ushort));
	for (ushort i = 0; i < nSpecies; ++i) fixedSpecArr[i] = i;

	//Fetch chemical env geoms (and set cartesian coordinates for imgs)
	printf("LOG: fetching chemical environment geometries ... ");
	for(ushort i = 0; i < nSites; ++i) MoveToCell(&allSites[i]);
	ushort imgArr[3]; GiveNImgs(&imgArr, rCut, A);
	env** envBase;
	SetSiteGeoms(nSites, &envBase, nSites, &allSites, imgArr, rCut, A, 
				 nSpecies);
	printf("ok\n");


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
	printf("ok\n");

	//Show absolute probabilities of sublattice species swaps
	//also overwrite default swapping probabilities if necessary
	for(ushort i = 0; i < nSublats; ++i){
		if(ovwSwapProbs) sublats[i]->swapProb = ovwSwapProbs[i];
		printf("JOB: p(s%u) = %f\n", i, sublats[i]->swapProb);
	}
	free(ovwSwapProbs);

	//Mapping of probabilities for uniform random choices of sublattices
	double* mapLo; double* mapHi;
	double* absPs = malloc(nSublats*sizeof(double));
	for(ushort i = 0; i < nSublats; ++i) absPs[i] = sublats[i]->swapProb;
	UniformMap(&mapLo, &mapHi, nSublats, absPs);
	free(absPs);

	//Avoid possible cataclysmic disaster
	if(N_REALLOC < 2*nSites + 1){
		printf("ERR: N_REALLOC < 2(number of sites) + 1\n");
		exit(1);
	}

	//Setup search tables, trees
	hbkt unqEnvTable[HASH_TABLE_LEN]; InitTable_h(&unqEnvTable);
	ushort envColls = 0;
	gbkt unqDecTable[GHASH_TABLE_LEN]; InitTable_gh(&unqDecTable);
	uint decColls = 0;


	//Structure generation-----------------------------------------------------
	//If a random site config vector results in a unique chem env vector, keep 
	//both of them

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
		uint nReppedEnvs = 0;
		ReadWarmStartFile(INFILE_WARM, &nUnqEnvC, &nUnqEnvM, &unqEnvs,
						  &nUnqDecC, &nUnqDecM, &unqDecs, &unqSites,
						  &nReppedEnvs, nSpecies, fixedSpecArr, nSites);
		nUnreppedEnvs = nDecsPerEnv*nUnqEnvC - nReppedEnvs;		

		////re-add things to search trees
		env* tmpEnv; uint dummy1 = 0; ushort dummy2 = 0; ushort dummy3 = 0;
		for(uint i = 0; i < nUnqEnvC; ++i){
			unqEnvs[i]->decId = i;
			tmpEnv = AddEnv(&unqEnvTable, nSpecies, unqEnvs[i], &envColls, 
							&dummy1);
		}
		for(uint i = 0; i < nUnqDecC; ++i){
			Add_gh(&unqDecTable, unqDecs[i], nSites, &decColls,
				   &dummy2, &dummy1, &dummy3);
		}

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
	////Initialize the current decomp array
	ushort* currDecArr = calloc(nUnqEnvM, sizeof(ushort));


	//Set all initial envs as needing to be recalced, do the first step
	for(ushort i = 0; i < nSites; ++i) envBase[i]->recalc = (ushort)1;

	///Continue to generate and eval random site permutations
	printf("GEN: n(tot)      n(unqS)  n(unqE)  n(lunqE)    n(unrE)"
		   "  n(collE) n(collD)\n");
	ushort swapSublat = 0; ushort swapSpecInd1 = 0; ushort swapSpecInd2 = 0;
	uvlong slue = 0llu; ///since last unique env
	env* oldEnv; env* newEnv; env* tmpE; uint* query;
	ushort hKey; uint addInd; ushort dSuccess;
	for(uvlong n = 0llu; n < nSteps; ++n){
		////set the current cell's decomposition...
		for(ushort i = 0; i < nSites; ++i){
			if(envBase[i]->recalc){
				envBase[i]->recalc = (ushort)0;

				/////Old env: the environment before resetting the comparison
				/////arrays - use to update the previous env decomposition
				/////On warm start, this won't mess anything up since 
				/////envBase[i]->sortedDists is set to garbage values for n = 0
				/////i.e. SearchEnv() will "always" return NULL if n = 0.  
				oldEnv = SearchEnv(unqEnvTable, envBase[i], nSpecies);
				if(oldEnv) currDecArr[oldEnv->decId]--; 
				SetSpecArr_e(&envBase[i], nSpecies);    
				SetDistArr_e(&envBase[i]);             
				/////New env: the updated environment, now consistent with 
				/////the random swaps - use to update the current env 
				/////decomposition, or add new envs
				newEnv = SearchEnv(unqEnvTable, envBase[i], nSpecies);
				
				/////If this new env already exists, just update the decomp
				/////counter
				if(newEnv) currDecArr[newEnv->decId]++;
				/////Otherwise its brand new - make a copy of it and store
				else{
					slue = 0llu;
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
					newEnv->sortedDists = malloc((newEnv->nDists)*
												  sizeof(double));
					SetDistArr_e(&newEnv);

					newEnv->decId = nUnqEnvC;
					newEnv->repped = (ushort)0;
					nUnreppedEnvs += nDecsPerEnv;

					unqEnvs[nUnqEnvC] = newEnv;
					tmpE = AddEnv(&unqEnvTable, nSpecies, newEnv, 
								  &envColls, &nUnqEnvC);
				}
			}
		}

		if(nUnreppedEnvs){
			////store sparse version of decomp + update repped envs if unique
			////and important
			unqDecs[nUnqDecC] = malloc(2*nSites*sizeof(uint));
			DecompConv(nSites, currDecArr, &unqDecs[nUnqDecC]);
			Add_gh(&unqDecTable, unqDecs[nUnqDecC], nSites, &decColls, 
				   &hKey, &addInd, &dSuccess);

			if(dSuccess){
				/////Update an env's number of unique decomp representations
				ushort accum = 0;
				for(ushort i = 0; accum < nSites; i += (ushort)2){
					if(unqEnvs[unqDecs[nUnqDecC][i]]->repped < nDecsPerEnv){
						unqEnvs[unqDecs[nUnqDecC][i]]->repped += (ushort)1;
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
				unqDecTable[hKey].nEntries--;
				unqDecTable[hKey].sparseArrs = realloc(
											   unqDecTable[hKey].sparseArrs, 
									unqDecTable[hKey].nEntries*sizeof(uint*)
													  );
				free(unqDecs[nUnqDecC]);
				if(unqDecTable[hKey].nEntries) decColls--;
			}
			else free(unqDecs[nUnqDecC]);
			
		}
		Keep: NOP

		////allocate more space if necessary
		if(nUnqEnvM - nUnqEnvC <= nSites){ /////wont overflow: max >= current
			nUnqEnvM += N_REALLOC;
			unqEnvs = realloc(unqEnvs, nUnqEnvM*sizeof(env*));
			currDecArr = realloc(currDecArr, nUnqEnvM*sizeof(ushort));
			/////need to make sure the extra decompositions are set to 0
			memset(currDecArr + nUnqEnvM - N_REALLOC, (ushort)0, 
				   sizeof(ushort)*N_REALLOC);
		}
		if(nUnqDecC == nUnqDecM){
			nUnqDecM += N_REALLOC;
			unqDecs = realloc(unqDecs, nUnqDecM*sizeof(uint*));
			unqSites = realloc(unqSites, nUnqDecM*sizeof(ushort*));
		}
		
		////make new random permutation, print out info
		UniformRandSwaps(nSublats, &sublats, mapLo, mapHi, 1u,
						 &swapSublat, &swapSpecInd1, &swapSpecInd2);
		MarkEnvsToRecalc(sublats, swapSublat, swapSpecInd1, swapSpecInd2);

#if GEN_PRINTSTYLE
		if(n%GEN_PRINTEVERY == 0llu){
#if (GEN_PRINTSTYLE == 1)
			printf("GEN: %011llu %08u %08u %011llu %08u %08u %08u\r", n, 
				   nUnqDecC, nUnqEnvC, slue, nUnreppedEnvs, envColls, 
				   decColls);
			fflush(stdout);
#endif
#if (GEN_PRINTSTYLE == 2)
			printf("GEN: %07.3f%% completed\r", (float)n/(float)nSteps*100.0f);
			fflush(stdout);
#endif
		}
#endif
		slue++;
	}
	printf("GEN: %011llu %08u %08u %011llu %08u %08u %08u\n", nSteps, nUnqDecC,
		   nUnqEnvC, slue, nUnreppedEnvs, envColls, decColls);

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
		free(allSites[i]->parentEnvs); ///free the actual envs later
		free(allSites[i]);
		for(ushort j = 0; j < envBase[i]->nSites; ++j){
			free(envBase[i]->sites[j]);
		}
		free(envBase[i]->sites);
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
		free(unqEnvs[i]->sites);
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
	DeallocTable_h(&unqEnvTable);
	DeallocTable_gh(&unqDecTable);

	return 0;
}