#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include "Decs.h"
#include "Constants.h"

#ifdef __unix__
#include "extern/kdtree-master/kdtree.h"
#endif
#ifdef _WIN32
#include "extern//kdtree-master//kdtree.h"
#endif

//Changes direct coordinates of a site s.t. 0 <= site coord i < 1
void MoveToCell(site** s){
	for(ushort i = 0; i < 3; ++i){
		while((*s)->crdsD[i] < 0.0) {
			(*s)->crdsD[i] += 1.0;
		}
		while((*s)->crdsD[i] >= 1.0) {
			(*s)->crdsD[i] -= 1.0;
		}
	}
	return;
}

//Returns the number of images necessary for periodic BCs.  
//We need enough periodic images s.t. rCut at any edge does not pass
//a periodic wall.  
void GiveNImgs(ushort nImgs[3], const double rCut, const double A[3][3]){
	for(unsigned int i = 0; i < 3; ++i){
		double latLen = sqrt(A[i][0]*A[i][0] +
							 A[i][1]*A[i][1] +
							 A[i][2]*A[i][2]);
		for(ushort j = 1;; ++j){
			if(rCut < (double)j*latLen){
				nImgs[i] = j;
				break;
			}
		}
	}
}

//Setup enviornment geometries taking periodic BCs into account
//Also attach parent envs to each important site
void SetSiteGeoms(const ushort nEnvs, env*** envs,
				  const ushort nSites, site*** sites,
				  const ushort nImgs[3], const double rCut, 
				  const double A[3][3], const ushort nSpecTot){

	///Get array of all neighbors 
	uint addSize = (2*nImgs[0]+1)*(2*nImgs[1]+1)*(2*nImgs[2]+1)*nSites;
	site** additionals = malloc(addSize*sizeof(site*));
	///Loop over all relevant images in neighboring cells
	uint count = 0;
	for(short a = -((short)nImgs[0]); a <= (short)nImgs[0]; ++a){
	for(short b = -((short)nImgs[1]); b <= (short)nImgs[1]; ++b){
	for(short c = -((short)nImgs[2]); c <= (short)nImgs[2]; ++c){
		for(ushort i = 0; i < nSites; ++i){
			site* thisSite = malloc(sizeof(site));
			for(ushort j = 0; j < 3; ++j){
				thisSite->crdsD[j] = (*sites)[i]->crdsD[j];
			}
			thisSite->species = (*sites)[i]->species;
			thisSite->self = (*sites)[i]->self;

			thisSite->crdsD[0] += (double)a;
			thisSite->crdsD[1] += (double)b;
			thisSite->crdsD[2] += (double)c;
			SetCarCrds_s(&thisSite, A);

			///-1: original, -2: periodic image
			///technically an overflow, yes, but who cares?  Nobody will have
			///this many elements in their cell anyways
			thisSite->species = (ushort)-2;
			if(a == 0u && b == 0u && c == 0u){
				thisSite->species = (ushort)-1;
			}

			additionals[count] = thisSite;
			count++;
		}
	}
	}
	}

	//Finish allocating
	for(ushort i = 0; i < nSites; ++i){
		(*sites)[i]->nParentEnvs = (ushort)0;
		(*sites)[i]->parentEnvs = malloc(0); ///NULL for now
	}

	//Connect sites to periodic images with a KD tree
	void* kdNode = kd_create(3);
	for(uint i = 0; i < addSize; ++i){
		kd_insert3(kdNode, additionals[i]->crdsC[0], additionals[i]->crdsC[1],
				   additionals[i]->crdsC[2], i);
	}

	///Set the list of enviornments
	*envs = malloc(nEnvs*sizeof(env*));
	double x = 0.0; double y = 0.0; double z = 0.0;
	ushort envCount = 0;
	for(uint i = 0; i < addSize; ++i){
		////Assign nbrs only to sites that lie within the orig cell
		if(additionals[i]->species != (ushort)-1){
			continue;
		}
		////temp variables - kd node CHANGES THE COORDINATES
		////somewhere in their functions...
		x = additionals[i]->crdsC[0];
		y = additionals[i]->crdsC[1];
		z = additionals[i]->crdsC[2];

		struct kdres* res = kd_nearest_range3(kdNode, x, y, z, rCut);
		int resSize = kd_res_size(res);

		env* thisEnv = malloc(sizeof(env));
		thisEnv->nSites = resSize;
		thisEnv->sites = malloc(resSize*sizeof(site*));
		AllocCmpArrs_e(&thisEnv, nSpecTot);

		////Assign neighbors to this enviornment,
		////assign parent envs to each site at the same time
		for(ushort j = 0; j < resSize; ++j){
			uint thisInd = (uint*)kd_res_item3(res, &x, &y, &z);

			site* nbr = malloc(sizeof(site));
			for(ushort k = 0; k < 3; ++k){
				nbr->crdsD[k] = additionals[thisInd]->crdsD[k];
				nbr->crdsC[k] = additionals[thisInd]->crdsC[k];
			}
			nbr->self = additionals[thisInd]->self;
			nbr->species = nbr->self->species;

			nbr->self->nParentEnvs++;
			nbr->self->parentEnvs = realloc(nbr->self->parentEnvs, 
										    (nbr->self->nParentEnvs)*
											sizeof(env*)
										   );			
			nbr->self->parentEnvs[nbr->self->nParentEnvs-(ushort)1] = thisEnv;
			

			thisEnv->sites[j] = nbr;
			kd_res_next(res);
		}


		(*envs)[envCount] = thisEnv;

		kd_res_free(res);
		envCount++;
	}

	kd_free(kdNode);
	for(uint i = 0; i < addSize; ++i){
		free(additionals[i]);
	}
	free(additionals);

	return;
}

//Sets a lattice to a specific configuration
void SetConfig(const ushort nSites, site*** sites, const ushort* species){
	for(ushort i = 0; i < nSites; ++i){
		*(*sites)[i]->species = species[i];
	}

	return;
}

//Takes a full array of environments (arr[a] = b means there are b instances
//of env # arr[a]) and creates a sparse vector (arr[a] = b, arr[a + 1] = c 
//means there are c instances of env # b ) for storage
void DecompConv(const ushort nTotEnvs, const ushort* full, uint **sparse){
	ushort accum = 0;
	ushort sparseLoc = 0;
	for(uint i = 0; accum < nTotEnvs; ++i){
		if(full[i]){
			(*sparse)[sparseLoc] = i;
			(*sparse)[sparseLoc + (ushort)1] = (uint)full[i];
			accum += full[i];
			sparseLoc += (ushort)2;
		}
	}

	return;
}

//Compares a full decomp vector to a sparse one
//Sparse vectors guarenteed to be written in ascending order wrt the full
//vector, so we don't need a double loop
int DecompCmpr_fs(const ushort nTotEnvs, 
				  const ushort* full, const uint* sparse){
	ushort accum = 0;
	for(uint itrS = 0; accum < nTotEnvs; itrS += (ushort)2){
		if(full[sparse[itrS]] != sparse[itrS + (ushort)1]) return 0;
		accum += sparse[itrS + (ushort)1];
	}

	return 1;
}

int DecompCmpr_ss(const ushort nTotEnvs, const uint* a, const uint* b){
	ushort accum = 0; ///only need one accum - if a(i+1] != b[i+1], we'll
					  ///return 0 anyways
	for(ushort itrS = 0; accum < nTotEnvs; itrS += (ushort)2){
		if(a[itrS] != b[itrS]) return 0;
		if(a[itrS + (ushort)1] != b[itrS + (ushort)1]) return 0;

		accum += (ushort)a[itrS + (ushort)1];
	}

	return 1;
}

//Deepcopies the current system state (the set of sublats, envs, sites)
//into a brand new state
//Requires a deepcopy of the fixed species array to have already been made
//Allocates all required memory
// *** NOT NECESSARY - JUST READ INPUT FILES TWICE
// *** KEEPING THIS HERE IN CASE I NEED IT LATER
// *** THE LOGIC IN THE LINE 'nbr->self = envsO[i]->sites[j]->self;' DOES NOT
//     WORK - IT MAKES A SHALLOW COPY
void DeepcopyState_DontUse(const ushort nSites, const site** sitesO, site*** sitesN,
				   const ushort nEnvs, const env** envsO, env*** envsN,
				   const ushort nSublats, 
				   const sublat** sublatsO, sublat*** sublatsN,
				   const ushort nSpecTot, const ushort* fixedSpecArr){
	
	//site-sublattice loop:
	//make new sites w/o parent envs, sublattices
	*sitesN = malloc(nSites*sizeof(site*));
	*sublatsN = malloc(nSublats*sizeof(sublat*));
	for(ushort i = 0, siteInd = 0; i < nSublats; ++i){
		(*sublatsN)[i] = malloc(sizeof(sublat));

		(*sublatsN)[i]->nSites = sublatsO[i]->nSites;
		(*sublatsN)[i]->swapProb = sublatsO[i]->swapProb;
		(*sublatsN)[i]->sites = malloc((sublatsO[i]->nSites)*sizeof(site*));
		for(ushort j = 0; j < sublatsO[i]->nSites; ++j, ++siteInd){

			///deepcopy site
			site* newSite = malloc(sizeof(site));
			newSite->species = malloc(sizeof(ushort));
			*newSite->species = *(sublatsO[i]->sites[j]->species);
			for(ushort k = 0; k < 3; ++k){
				newSite->crdsC[k] = sublatsO[i]->sites[j]->crdsC[k];
				newSite->crdsD[k] = sublatsO[i]->sites[j]->crdsD[k];
			}
			newSite->self = newSite; 

			////prep for setting parent envs
			newSite->nParentEnvs = (ushort)0;
			newSite->parentEnvs = malloc(0);

			(*sitesN)[siteInd] = newSite;
			(*sublatsN)[i]->sites[j] = newSite;
		}
	}

	//site-env loop:
	//make new envs, use previously made sites, add parent envs
	*envsN = malloc(nEnvs*sizeof(env*));
	for(ushort i = 0; i < nEnvs; ++i){

		///deepcopy env
		env* newEnv = malloc(sizeof(env));
		newEnv->nSites = envsO[i]->nSites;
		newEnv->sites = malloc(newEnv->nSites*sizeof(site*));
		AllocCmpArrs_e(&newEnv, nSpecTot);

		////new neighbor
		for(ushort j = 0; j < newEnv->nSites; ++j){
			site* nbr = malloc(sizeof(site));
			for(ushort k = 0; k < 3; ++k){
				nbr->crdsD[k] = envsO[i]->sites[j]->crdsD[k];
				nbr->crdsC[k] = envsO[i]->sites[j]->crdsC[k];
			}
			nbr->self = envsO[i]->sites[j]->self;
			nbr->species = nbr->self->species;

			nbr->self->nParentEnvs++;
			nbr->self->parentEnvs = realloc(nbr->self->parentEnvs,
											(nbr->self->nParentEnvs)*
											sizeof(env*)
										   );
			nbr->self->parentEnvs[nbr->self->nParentEnvs-(ushort)1] = newEnv;

			newEnv->sites[j] = nbr;
		}

		(*envsN)[i] = newEnv;
	}
}

void DeepcopyState(const ushort nSites, const site** src, site*** dst){
	for(ushort i = 0; i < nSites; ++i){
		*(*dst)[i]->species = *src[i]->species;
	}
}