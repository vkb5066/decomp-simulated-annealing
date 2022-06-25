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
void SetSiteGeoms(const ushort nEnvs, env*** envs,
				  const ushort nSites, const site** sites,
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
				thisSite->crdsD[j] = sites[i]->crdsD[j];
			}
			thisSite->species = sites[i]->species;
			thisSite->self = sites[i]->self;

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

		////Assign neighbors to this enviornment
		for(ushort j = 0; j < resSize; ++j){
			uint thisInd = (uint*)kd_res_item3(res, &x, &y, &z);

			site* nbr = malloc(sizeof(site));
			for(ushort k = 0; k < 3; ++k){
				nbr->crdsD[k] = additionals[thisInd]->crdsD[k];
				nbr->crdsC[k] = additionals[thisInd]->crdsC[k];
			}
			nbr->self = additionals[thisInd]->self;
			nbr->species = nbr->self->species;
			

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
//Sparse vectors guarenteed to be written in ascendnig order wrt the full
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