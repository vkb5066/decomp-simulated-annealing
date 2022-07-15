#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "Decs.h"
#include "Constants.h"

//Site Functions---------------------------------------------------------------
//Set Cartesian coordinates given direct coordinates and the lattice
void SetCarCrds_s(site** s, const double A[3][3]){
	for(ushort i = 0; i < 3; ++i){
		(*s)->crdsC[i] = 0.0;
		for(ushort j = 0; j < 3; ++j){
			(*s)->crdsC[i] += (*s)->crdsD[j] * A[j][i];
		}
	}
	return;
}

//Squared distance between two sites
double DSqd_s(const site* a, const site* b) {
	return(
		(b->crdsC[0] - a->crdsC[0]) * (b->crdsC[0] - a->crdsC[0]) +
		(b->crdsC[1] - a->crdsC[1]) * (b->crdsC[1] - a->crdsC[1]) +
		(b->crdsC[2] - a->crdsC[2]) * (b->crdsC[2] - a->crdsC[2]) +
		((double)*b->species - (double)*a->species) *
		((double)*b->species - (double)*a->species)
		  );
}

//Prints species followed by direct coordinates then cartesian coordinates
void Print_s(const site* s){
	printf("%u %.11f %.11f %.11f %.11f %.11f %.11f", *s->species, 
		   s->crdsD[0], s->crdsD[1], s->crdsD[2],
		   s->crdsC[0], s->crdsC[1], s->crdsC[2]);
	return;
}


//Env Functions----------------------------------------------------------------
//Sorting function
int dblcmpr(const void* a, const void* b){
	if (*(const double*)a < *(const double*)b) return -1;
	if (*(const double*)a > *(const double*)b) return +1;
	return 0;
}

//Allocate space for species and distance arrays
//The geometry of an env does not change, so the distance array has a fixed
//size.  The species array must be reset to {0}, though, on every element swap
void AllocCmpArrs_e(env** e, const ushort nSpecTot){
	(*e)->sortedSpecies = malloc(nSpecTot*sizeof(ushort));
	(*e)->nDists = (ushort)((*e)->nSites)*((*e)->nSites - 1u)/2u;
	(*e)->sortedDists = malloc(((*e)->nDists)*sizeof(double));
}

//Sets an array where the ith component provides the count of how 
//many elements of type i are in this enviornment
void SetSpecArr_e(env** e, const ushort nSpecTot){
	memset((*e)->sortedSpecies, 0, nSpecTot*sizeof(ushort));
	for(ushort i = 0; i < (*e)->nSites; ++i){
		(*e)->sortedSpecies[*((*e)->sites[i]->species)] += (ushort)1;
	}
	return;
}

//Sets an array of unique inter-site distances, sorted
void SetDistArr_e(env** e){
	for(ushort i = 0u, loc = 0u; i < (*e)->nSites - 1; ++i){
		for(ushort j = i + 1u; j < (*e)->nSites; ++j, ++loc){
			(*e)->sortedDists[loc] = DSqd_s((*e)->sites[i], (*e)->sites[j]);
		}
	}
#if (SORT_TYPE == 0)
	qsort((*e)->sortedDists, (*e)->nDists, sizeof(double), dblcmpr);
#endif
#if (SORT_TYPE == 1)
	ISort_d(&((*e)->sortedDists), (*e)->nDists);
#endif
	return;
}

//Returns 1 if two environments are "equal", as measured by their 4D distance
//tables (see DSqrd_s)
int Cmpr_e(const env* e1, const env* e2, const ushort nSpecTot){
	//Most likley discrep: unequal # of elements
	for(ushort i = 0; i < nSpecTot; ++i){
		if(e1->sortedSpecies[i] != e2->sortedSpecies[i]) return 0;
	}

	//Second most likley discrep: unequal number of sites
	if(e1->nDists != e2->nDists) return 0;

	//Longest check: look over all distances table for discreps
	for(ushort i = 0; i < e1->nDists; ++i){
		if(fabs(e1->sortedDists[i] - e2->sortedDists[i]) > CART_EPS) return 0;
	}

	return 1;
}

//Swaps the data needed for comparing two envs i.e. the fixed species arrays
//and the sorted distance arrays
//Be sure that src and dst don't point to the same memory location!
void SwapCmpArrs_e(const env*restrict src, env*restrict*restrict dst, 
				   const ushort nSpecTot){
	for(ushort i = 0; i < nSpecTot; ++i){
		(*dst)->sortedSpecies[i] = src->sortedSpecies[i];
	}
	for(ushort i = 0; i < src->nDists; ++i){
		(*dst)->sortedDists[i] = src->sortedDists[i];
	}
	return;
}

//Prints an environment as nSites followed by [nSites] site prints
//if give* = 0, prints out NULL instead
void Print_e(const env* e, ushort giveNrg, ushort giveRep){
	printf("%u\n", e->nSites);
	for(ushort i = 0; i < e->nSites; ++i){
		Print_s(e->sites[i]); 
		printf("\n");
	}
	if(giveNrg) printf("%.11f\n", e->nrg);
	else printf("NULL\n");
	if(giveRep) printf("%u", e->repped);
	else printf("NULL");

	return;
}


//Env tree functions-----------------------------------------------------------
void Alloc_b(branch** b, double key){
	*b = malloc(sizeof(branch));
	(*b)->thisEnv = NULL;
	(*b)->key = key;
	(*b)->nBranches = (ushort)0;
	(*b)->branches = malloc(0);
}

void Dealloc_b(branch** b){
	for(ushort i = 0; i < (*b)->nBranches; ++i){
		Dealloc_b(&((*b)->branches[i]));
	}

	free((*b)->branches);
	free(*b); ///do NOT free the env pointers those are done elsewhere
}

//Adds a pointer to an env to a tree.  Increments (NOT SETS!) a success counter
//if a new env was added.  Returns a pointer to either (a) the env that was
//sucessfully added or (b) the env that conflicts with the env that was 
//attempted to be added
env* Add_b(branch*restrict*restrict head, const env*restrict toAdd, 
		   uint*restrict success){
	branch* itr = *head;
	for(ushort i = 0; i < toAdd->nDists; ++i){
		///check to see if we have a path to walk down
		for(ushort j = 0; j < itr->nBranches; ++j){
			////if so, go down it
			if(fabs(toAdd->sortedDists[i] - itr->branches[j]->key) < CART_EPS){
				itr = itr->branches[j];
				goto Walked;
			}
		}
		
		////otherwise, make the path and then go down it
		itr->nBranches++;
		itr->branches = realloc(itr->branches, itr->nBranches*sizeof(branch*));
		branch* thisBranch;
		Alloc_b(&thisBranch, toAdd->sortedDists[i]);
		itr->branches[itr->nBranches - (ushort)1] = thisBranch;
		itr = thisBranch;

		Walked: NOP
	}

	//Now, itr is at the end of the tree
	///Return the env pointer if an env exists here already
	if(itr->thisEnv) return itr->thisEnv;
	
	///If we're here, an env does not exist
	itr->thisEnv = toAdd;
	(*success)++;
	return itr->thisEnv;
}

//Searches for an env in the env tree.  Returns a pointer to the env 
//that was sucessfully found, or NULL if it was not found
env* Search_b(const branch*restrict head, const env*restrict toFind){
	branch* itr = head;
	for(ushort i = 0; i < toFind->nDists; ++i){
		///check to see if we have a path to walk down
		for(ushort j = 0; j < itr->nBranches; ++j){
			////if so, go down it
			if(fabs(toFind->sortedDists[i]-itr->branches[j]->key) < CART_EPS){
				itr = itr->branches[j];
				goto Walked;
			}
		}
		///if we're here, we can't find a path: return NULL
		return NULL;

		Walked: continue;
	}
	return itr->thisEnv;
}

//Hash Functions---------------------------------------------------------------

//Hashes an array to an int: as suggested by this question on SO
//.../2351087/what-is-the-best-32bit-hash-function-for-short-strings-tag-names
//and answered by Nick Dandoulakis.  
ushort ArrToHash_h(const ushort len, const ushort* arr){
	ushort bkt = arr[0];
	for(ushort i = 1; i < len; ++i){
		bkt = bkt*HASH_TABLE_MUL + arr[i];
	}
	return bkt%HASH_TABLE_LEN;
}

ushort ArrToHash_gh(const ushort len, const ushort* arr){
	ushort bkt = arr[0];
	for(ushort i = 0; i < len; ++i){
		bkt = bkt*GHASH_TABLE_MUL + arr[i];
	}
	return bkt%GHASH_TABLE_LEN;
}

//Initializes a hash table of size HASH_TABLE_LEN and returns a pointer
//to it
void InitTable_h(hbkt* ptr){
	for(ushort i = 0; i < HASH_TABLE_LEN; ++i){
		ptr[i].nEntries = (ushort)0;
		ptr[i].specArrs = malloc(0);
		ptr[i].branchHeads = malloc(0);
	}
	return;
}

void InitTable_gh(gbkt* ptr){
	for(ushort i = 0; i < GHASH_TABLE_LEN; ++i){
		ptr[i].nEntries = 0u;
		ptr[i].sparseArrs = malloc(0);
	}
	return;
}

void DeallocTable_h(hbkt* h){
	for(ushort i = 0; i < HASH_TABLE_LEN; ++i){
		for(ushort j = 0; j < h[i].nEntries; ++j){
			free(h[i].specArrs[j]);
			Dealloc_b(&(h[i].branchHeads[j]));
		}
		free(h[i].specArrs);
		free(h[i].branchHeads);
	}
}

void DeallocTable_gh(gbkt* g){
	for(ushort i = 0; i < GHASH_TABLE_LEN; ++i){
		free(g[i].sparseArrs); ///will free the actual sparse arrays elsewhere
	}
}

//Adds an entry to a hash table.  Counts increments a collision counter if any
//collisions happen
//Returns a pointer to (a) the branch that was added or (b) the branch that 
//exists in the spot where the specArr was meant to be added
branch* Add_h(hbkt*restrict table, const ushort hashInd,
			  const ushort specArrLen, const ushort*restrict specArr,
			  ushort*restrict collisions){
	hbkt* bkt = table + hashInd;

	//Check for identical entries
	for(ushort i = 0; i < bkt->nEntries; ++i){
		ushort arrsEqu = 1;
		for(ushort j = 0; j < specArrLen; ++j){
			if(specArr[j] != bkt->specArrs[i][j]){
				arrsEqu = (ushort)0;
				break;
			}
		}
		if(arrsEqu) return bkt->branchHeads[i];
	}

	//If we're here, we haven't found an identical entry - add a new one
	///Add species array
	bkt->specArrs = realloc(bkt->specArrs, 
						    (bkt->nEntries + (ushort)1)*sizeof(ushort*));
	bkt->specArrs[bkt->nEntries] = malloc(specArrLen*sizeof(ushort));
	for(ushort i = 0; i < specArrLen; ++i){
		bkt->specArrs[bkt->nEntries][i] = specArr[i];
	}
	///Add pointer to new tree head node
	bkt->branchHeads = realloc(bkt->branchHeads, 
							   (bkt->nEntries + (ushort)1)*sizeof(branch*));
	Alloc_b(&(bkt->branchHeads[bkt->nEntries]), 0.0);

	bkt->nEntries++;
	if(bkt->nEntries > (ushort)1) (*collisions)++;

	return bkt->branchHeads[bkt->nEntries - (ushort)1];	
}

//Adds a sparse array to this table index
//if the entry was sucessfully added, sets suc = 1 and fills the hash and
//addInd values with the hash key and the index of the added sparse arr
//if a duplicate entry was found, sets suc = 0
void Add_gh(gbkt*restrict table, const uint*restrict sparseArr,
			 const ushort nTotEnvs, uint*restrict colls,
			 ushort*restrict hash, uint*restrict addInd, ushort*restrict suc){
	ushort slen = 0; ushort accum = 0;
	while(accum < nTotEnvs){
		accum += sparseArr[slen + (ushort)1];
		slen += (ushort)2;
	}

	ushort hKey = ArrToHash_gh(slen, sparseArr);
	gbkt* bkt = table + hKey;

	///Check for identical entries
	for(uint i = 0; i < bkt->nEntries; ++i){
		if(DecompCmpr_ss(nTotEnvs, sparseArr, bkt->sparseArrs[i])){
			(*suc) = (ushort)0;
			return;
		}
	}

	///If we're here, we havent found any identical entries - set the new one
	bkt->sparseArrs = realloc(bkt->sparseArrs,
							  (bkt->nEntries + 1u)*sizeof(uint*));
	bkt->sparseArrs[bkt->nEntries] = sparseArr;
	(*hash) = hKey;
	(*addInd) = bkt->nEntries;
	bkt->nEntries++;
	if(bkt->nEntries > 1u) (*colls)++;
	(*suc) = (ushort)1;

	return;
}


//Searches for a hash table entry given an array of shorts
//Returns a pointer to the branch that corresponds to the search array
//and NULL if no such entry exists
branch* Search_h(const hbkt*restrict table, const ushort hashInd, 
				 const ushort specArrLen, const ushort*restrict specArr){
	hbkt* bkt = table + hashInd;

	if(!bkt->nEntries) return NULL;
#if !(SAFE_HASH_LOOKUP)
	if(bkt->nEntries == (ushort)1) return bkt->branchHeads[0];
#endif
	for(ushort i = 0; i < bkt->nEntries; ++i){
		ushort arrsEqu = 1;
		for(ushort j = 0; j < specArrLen; ++j){
			if(specArr[j] != bkt->specArrs[i][j]){
				arrsEqu = (ushort)0;
				break;
			}
		}
		if(arrsEqu) return bkt->branchHeads[i];
	}

	return NULL;
}


//Adds an env pointer to the hash-table-to-lookup-tree monstrosity
//Increments a collision counter (for the hash table) and a counter
//if the environment to add is unique
env* AddEnv(hbkt*restrict table, const ushort nSpecTot, 
			const env*restrict toAdd,
			ushort*restrict collCnt, uint*restrict newEnvCount){
	//Get the pointer to the lookup tree from hash table
	branch* thisBranch = Add_h(&(*table), 
							   ArrToHash_h(nSpecTot, toAdd->sortedSpecies), 
							   nSpecTot, toAdd->sortedSpecies, 
							   &(*collCnt));
	return Add_b(&thisBranch, toAdd, &(*newEnvCount));
}

//Searches the hash-table-to-lookup-tree monstrosity for the query env.
//If found, returns a pointer to the corresponding env.  Otherwise returns NULL
env* SearchEnv(const hbkt*restrict table, const env*restrict query, 
			   const ushort nSpecTot){
	branch* thisBranch = Search_h(table, 
								  ArrToHash_h(nSpecTot, query->sortedSpecies), 
								  nSpecTot, query->sortedSpecies);
	if(!thisBranch) return NULL;
	return Search_b(thisBranch, query);
}