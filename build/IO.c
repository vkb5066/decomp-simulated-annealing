#pragma warning(disable:4996)
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "Decs.h"
#include "Constants.h"

//Input------------------------------------------------------------------------

//Reads A, site mapping, and lattice from a file
void ReadLatticeFile(const char* fileName,
					 double A[3][3],
					 ushort* nSpecies,
					 ushort* nSites, site*** sites,
					 ushort* nSublats,
					 sublat*** sublats){

	FILE* infile;
	infile = fopen(fileName, "r");
	if(!infile){
		printf("\nERR: unable to open %s\n", fileName);
		exit(1);
	}


	char line[LINESIZE];
	char* next;
	ushort readLat = 0;
	ushort readSit = 0;
	for(uint rlCount = 0; (readLat + readSit) < 2u; ++rlCount){
		if(rlCount > INFILE_LINE_MAX){
			printf("\nERR: missing keyword(s) in %s (or file is too long)\n", 
				   fileName);
			fclose(infile);
			exit(1);
		}
		fgets(line, LINESIZE, infile);

		//Unit / supercell definition
		if(strstr(line, "BEGIN CELL")){
			for(ushort i = 0; i < 3; ++i){
				fgets(line, LINESIZE, infile);
				A[i][0] = strtod(line, &next);
				A[i][1] = strtod(next, &next);
				A[i][2] = strtod(next, NULL);
			}
			readLat = 1;
			continue;
		}

		//Sublattice definition (oh boy)
		if(strstr(line, "BEGIN SUBLATTICE")){
			///total number of species in cell
			fgets(line, LINESIZE, infile);
			*nSpecies = (ushort)strtol(line, &next, BASE);
			char** speciesArr = malloc((*nSpecies)*sizeof(char*));
			for(ushort i = 0; i < (*nSpecies); ++i){
				speciesArr[i] = malloc(ELEMSIZE*sizeof(char));
				char tmp[ELEMSIZE];
				fscanf(infile, "%s", &tmp);
				for(ushort j = 0; j < ELEMSIZE; ++j){
					speciesArr[i][j] = tmp[j];
					if(tmp[j] == '\0'){ 
						break;
					}
				}
			}

			///Sublattice info
			fgets(line, LINESIZE, infile); ////endline char
			fgets(line, LINESIZE, infile);
			*nSublats = (ushort)strtol(line, &next, BASE);
			*sublats = malloc((*nSublats)*sizeof(sublat*));
			(*nSites) = 0;
			ushort* specCount = calloc(*nSpecies, sizeof(ushort));
			///loop over all sublattices ...
			for(ushort i = 0; i < (*nSublats); ++i){
				////first read in species info for this sublattice
				fgets(line, LINESIZE, infile);
				ushort nSpecs = (ushort)strtol(line, &next, BASE);
				ushort nSites_ = 0;
				for(ushort j = 0; j < nSpecs; ++j){
					char spec[ELEMSIZE]; int cnt;
					fscanf(infile, "%s", &spec);
					fscanf(infile, "%u", &cnt);
					for(ushort k = 0; k < *nSpecies; ++k){
						if(!strcmp(spec, speciesArr[k])){
							specCount[k] = (ushort)cnt;
							nSites_+= cnt;
							(*nSites) += cnt;
							break;
						}
					}
				}

				////Make new sublattice
				sublat* thisSublat = malloc(sizeof(sublat));
				thisSublat->nSites = nSites_;
				thisSublat->sites = malloc(nSites_*sizeof(site*));
				(*sublats)[i] = thisSublat;

				////set unscaled probabilities
				///-1.0 so that pswap(sublats w/ only one site, spec) = 0 
				(*sublats)[i]->swapProb = CountCombos(*nSpecies, 
													  specCount, 
													  nSites_) - 1.0;

				////then read in all of the coordinates for this sublattice
				fgets(line, LINESIZE, infile); ////endline char
				for(ushort j = 0; j < nSites_; ++j){
					fgets(line, LINESIZE, infile);

					/////Make new site
					site* thisSite = malloc(sizeof(site));
					thisSite->species = malloc(sizeof(ushort));
					for(ushort k = 0; k < *nSpecies; ++k){
						if(specCount[k]){
							*thisSite->species = k;
							specCount[k]--;
							break;
						}
					}
					thisSite->crdsD[0] = strtod(line, &next);
					thisSite->crdsD[1] = strtod(next, &next);
					thisSite->crdsD[2] = strtod(next, NULL);
					thisSite->self = thisSite;

					(*sublats)[i]->sites[j] = thisSite;
				}
			}

			///Normalize swapping probabilities
			double norm = 0.0;
			for(ushort i = 0; i < *nSublats; ++i){
				norm += (*sublats)[i]->swapProb;
			}
			for(ushort i = 0; i < *nSublats; ++i){
				(*sublats)[i]->swapProb /= norm;
			}

			for(ushort i = 0; i < *nSpecies; ++i){
				free(speciesArr[i]);
			}
			free(speciesArr);
			free(specCount);
			readSit = 1;
			continue;
		}

	}
	
	//Fill an array of sites that are not seperated into sublattices for easy
	//freeing, writing, etc.
	(*sites) = malloc(*nSites*sizeof(site*));
	for(ushort i = 0, cnt = 0; i < *nSublats; ++i){
		for(ushort j = 0; j < (*sublats)[i]->nSites; ++j, ++cnt){
			(*sites)[cnt] = (*sublats)[i]->sites[j];
		}
	}

	fclose(infile);
	return;
}

//Reads job parameters i.e. cutoff rad, s.a. stuff, ...
ushort ReadParamFile(const char* fileName, double* rCut,
					 uvlong* nSteps, uvlong* startStep,
					 double* startNrg, int* rSeed, 
					 ushort** initState){

	FILE* infile;
	infile = fopen(fileName, "r");
	if(!infile){
		printf("\nERR: Unable to open %s\n", fileName);
		exit(1);
	}


	char line[LINESIZE];
	char* next;
	ushort readRad = 0;
	ushort readAnn = 0;
	ushort readCon = 0;
	for(uint rlCount = 0; (readRad + readAnn + readCon) < 3u; ++rlCount){
		///the initial configuration is optional
		if(rlCount > INFILE_LINE_MAX && readRad && readAnn){
			break;
		}
		if(rlCount > INFILE_LINE_MAX){
			printf("\nERR: missing keyword(s) in %s (or file is too long)\n",
				   fileName);
			fclose(infile);
			exit(1);
		}
		fgets(line, LINESIZE, infile);

		//Cutoff radius
		if(strstr(line, "BEGIN CUTOFF")){
			fgets(line, LINESIZE, infile);
			*rCut = strtod(line, &next);
			
			readRad = 1;
			continue;
		}

		//Annealing params
		if(strstr(line, "BEGIN ANNEAL")){
			fgets(line, LINESIZE, infile);
			*nSteps = strtoull(line, &next, BASE); 
			*startStep = strtoull(next, NULL, BASE);
			fgets(line, LINESIZE, infile);
			*startNrg = strtod(line, NULL);
			fgets(line, LINESIZE, infile);
			*rSeed = strtol(line, &next, BASE);
			
			readAnn = 1;
			continue;
		}

		//Initial configuration
		if(strstr(line, "BEGIN STARTCONFIG")){
			fgets(line, LINESIZE, infile);
			ushort nSites = (ushort)strtol(line, &next, BASE);
			//fgets(line, LINESIZE, infile); ///move to next line
			for(ushort i = 0; i < nSites; ++i){
				int tmp = -1;
				fscanf(infile, "%u", &tmp);
				(*initState)[i] = (ushort)tmp;
			}

			readCon = 1;
		}

	}

	fclose(infile);
	return readCon;
}


void ReadEnvFile(const char* fileName, 
				 uint* nEnvs, env*** envs,
				 const ushort nSpecTot, const ushort* fixSpecArr){

	FILE* infile;
	infile = fopen(fileName, "r");
	if(!infile){
		printf("\nERR: unable to open %s\n", fileName);
		exit(1);
	}


	char line[LINESIZE];
	char* next;
	ushort readEnv = 0;
	ushort readDec = 0;
	for(uint rlCount = 0;; ++rlCount){
		if(rlCount > INFILE_LINE_MAX && readEnv){
			break;
		}
		if(rlCount > INFILE_LINE_MAX){
			printf("\nERR: missing keyword(s) in %s (or file is too long)\n",
				   fileName);
			fclose(infile);
			exit(1);
		}
		fgets(line, LINESIZE, infile);

		//Read chemical environments
		if(strstr(line, "BEGIN ENVS")){
			fgets(line, LINESIZE, infile);
			*nEnvs = (uint)strtoul(line, &next, BASE);
			*envs = malloc((*nEnvs)*sizeof(env*));
			for(uint i = 0; i < *nEnvs; ++i){
				env* thisEnv = malloc(sizeof(env));

				fgets(line, LINESIZE, infile);
				thisEnv->nSites = (ushort)strtoul(line, &next, BASE);
				thisEnv->sites = malloc((thisEnv->nSites)*sizeof(site*));
				AllocCmpArrs_e(&thisEnv, nSpecTot);
				for(ushort j = 0; j < thisEnv->nSites; ++j){
					fgets(line, LINESIZE, infile);
					site* thisSite = malloc(sizeof(site));
					
					ushort spec = (ushort)strtoul(line, &next, BASE);
					thisSite->species = fixSpecArr + spec; ///<- points to cnst
					thisSite->crdsD[0] = strtod(next, &next); ///species array:
					thisSite->crdsD[1] = strtod(next, &next); ///no need to 
					thisSite->crdsD[2] = strtod(next, &next); ///allocate a ton
					thisSite->crdsC[0] = strtod(next, &next); ///of extra ints 
					thisSite->crdsC[1] = strtod(next, &next); ///if theyre  
					thisSite->crdsC[2] = strtod(next, NULL); ///fixed
		
					thisEnv->sites[j] = thisSite;          
				}
				SetSpecArr_e(&thisEnv, nSpecTot);
				SetDistArr_e(&thisEnv);

				fgets(line, LINESIZE, infile);
				thisEnv->nrg = strtod(line, NULL);
				fgets(line, LINESIZE, infile);
				thisEnv->repped = strtoul(line, NULL, BASE);
				(*envs)[i] = thisEnv;
			}

			readEnv = 1;
		}

	}

	fclose(infile);
	return;
}

//Reads warm start file for environments.
//Fixed format: First envs, then env decompositions, then site writeouts
//envs: line1 = nEnvs, followed by those envs
//decs: line1 = number of rows, followed by those rows (sparse format i.e.
//cvcvcv ... : c = col = the env index
//             v = val = the number of envs of left column
//each row has at max, nSites c's - so 2*nSites entries.
//The actual size of the row is determined by the v's: when their sum (moving 
//leftward) = nSites, we must have run out of representations
//sites: line1 = nRows, followed by those rows (dense format of constant width
//nSites)
void ReadWarmStartFile(const char* fileName, 
					   uint* cSizeE, uint* mSizeE, env*** arrE,
					   uint* cSizeD, uint* mSizeD, uint*** arrD,
					   ushort*** arrS,
					   uint* nReppedEnvs,
					   const ushort nSpecTot, const ushort* fixSpecArr,
					   const ushort nSites){
	FILE* infile;
	infile = fopen(fileName, "r");
	if (!infile) {
		printf("\nERR: unable to open %s\n", fileName);
		exit(1);
	}

	char line[LINESIZE];
	char* next;

	//Fill array with envs
	fgets(line, LINESIZE, infile);
	*cSizeE = (uint)strtoul(line, NULL, BASE);
	*mSizeE = 0u;
	while((int)*mSizeE - (int)*cSizeE <= (int)nSites) *mSizeE += N_REALLOC;
	*arrE = malloc((*mSizeE)*sizeof(env*));
	*nReppedEnvs = 0;
	for(uint i = 0; i < *cSizeE; ++i){
		env* thisEnv = malloc(sizeof(env));

		fgets(line, LINESIZE, infile);
		thisEnv->nSites = (ushort)strtoul(line, &next, BASE);
		thisEnv->sites = malloc((thisEnv->nSites)*sizeof(site*));
		AllocCmpArrs_e(&thisEnv, nSpecTot);
		for(ushort j = 0; j < thisEnv->nSites; ++j){
			fgets(line, LINESIZE, infile);
			site* thisSite = malloc(sizeof(site));

			ushort spec = (ushort)strtoul(line, &next, BASE);
			thisSite->species = fixSpecArr + spec; ///<- points to cnst
			thisSite->crdsD[0] = strtod(next, &next); ///species array:
			thisSite->crdsD[1] = strtod(next, &next); ///no need to 
			thisSite->crdsD[2] = strtod(next, &next); ///allocate a ton
			thisSite->crdsC[0] = strtod(next, &next); ///of extra ints 
			thisSite->crdsC[1] = strtod(next, &next); ///if theyre  
			thisSite->crdsC[2] = strtod(next, NULL); ///fixed

			thisEnv->sites[j] = thisSite;
		}
		SetSpecArr_e(&thisEnv, nSpecTot);
		SetDistArr_e(&thisEnv);

		fgets(line, LINESIZE, infile); ///energy: skip - don't need
		fgets(line, LINESIZE, infile);
		thisEnv->repped = strtoul(line, NULL, BASE);
		(*nReppedEnvs) += thisEnv->repped;
		(*arrE)[i] = thisEnv;
	}


	//Fill array with env decompositions
	fgets(line, LINESIZE, infile);
	*cSizeD = (uint)strtoul(line, NULL, BASE);
	*mSizeD = 0u;
	while(*mSizeD <= *cSizeD) *mSizeD += N_REALLOC;
	*arrD = malloc((*mSizeD)*sizeof(uint*));
	for(uint i = 0; i < *cSizeD; ++i){
		(*arrD)[i] = malloc(2*nSites*sizeof(uint));
		ushort accum = 0;
		for(ushort j = 0; accum < nSites; j += (ushort)2){
			uint c = -1; uint v = -1;
			fscanf(infile, "%u", &c); ///env ID number
			(*arrD)[i][j] = c;
			fscanf(infile, "%u", &v);
			(*arrD)[i][j + (ushort)1] = v;

			accum += (ushort)((*arrD)[i][j + (ushort)1]);
		}

	}
	fgets(line, LINESIZE, infile); ///skip newline char


	//Fill array with site lists
	fgets(line, LINESIZE, infile);	  ///dummy - we already known nRows from 
	*arrS = malloc((*mSizeD)*sizeof(ushort*)); ///the env decomp row numbers
	for(uint i = 0; i < *cSizeD; ++i){
		(*arrS)[i] = malloc(nSites*sizeof(ushort));
		for(ushort j = 0; j < nSites; ++j){
			uint spec = -1;
			fscanf(infile, "%u", &spec);
			(*arrS)[i][j] = (ushort)spec;
		}
	}
	
	fclose(infile);
	return;
}


//Output-----------------------------------------------------------------------

void SitePrintout(const uint nPrints, const ushort nSites, const ushort** arr){
	printf("%u %u\n", nPrints, nSites);
	for(uint i = 0; i < nPrints; ++i){
		for(ushort j = 0; j < nSites; ++j){
			printf("%u ", arr[i][j]);
		}
		printf("\n");
	}
}

void EnergyPrintout(const uint nPrints, const double* arr){
	printf("%u\n", nPrints);
	for(uint i = 0; i < nPrints; ++i){
		printf("%f\n", arr[i]);
	}
}

void DecPrintout(const uint nPrints, const uint maxCol, 
				 const ushort nSites, const uint** arr){
	printf("%u %u\n", nPrints, maxCol);
	for(uint i = 0; i < nPrints; ++i){
		ushort accum = 0;
		for(ushort itrS = 0; accum < nSites; itrS += (ushort)2){
			printf("%u %u ", arr[i][itrS], arr[i][itrS + (ushort)1]);
			accum += (ushort)arr[i][itrS + (ushort)1];
		}
		printf("\n");
	}
}

void EnvPrintout(const uint nEnvs, const env** envs,
				 ushort giveNrg, ushort giveRep){
	printf("%u\n", nEnvs);
	for(uint i = 0; i < nEnvs; ++i) {
		Print_e(envs[i], giveNrg, giveRep);
		printf("\n");
	}
}

//Handles situations where we can't find an environment with a given energy
//Prints out the configuration, step number, number of unknowns, and 
//unknown environments in question
void UnknownEnvHandler(const ushort nSites, const site** sites,
					   const uvlong step, 
					   const ushort nUnknowns, const ushort* unknownEnvs,
					   const ushort nTrialEnvs, const env** trialEnvs){
	printf("ERR: at step %llu, the following config had %u missing envs:\n", 
		   step, nUnknowns);
	printf("ERR: BEGIN MISSING SITE PRINTOUT\n");
	printf("%u %u\n", 1u, nSites);
	for(ushort i = 0; i < nSites; ++i) printf("%u ", *(sites[i]->species));
	printf("\nEND missing SITE PRINTOUT\n");
	printf("ERR: BEGIN MISSING ENV PRINTOUT:\n");
	printf("%u\n", nUnknowns);
	for(ushort i = 0; i < nTrialEnvs; ++i){
		if(unknownEnvs[i]){
			Print_e(trialEnvs[i], (ushort)0, (ushort)0);
			printf("\n");
		}
	}
	printf("ERR: END MISSING ENV PRINTOUT");
	return;
}
