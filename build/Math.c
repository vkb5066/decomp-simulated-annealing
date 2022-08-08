#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "Decs.h"

//General math functions

//n choose k for large numbers (returns double since I'll be dividing w/ it)
double NCK(const uint n, const uint k){
	if(k == 0u || k == n){
		return 1.0;
	}
	uint k_ = (n-k < k)? n-k : k;
	
	double res = 1.0;
	for(uint i = 0u; i < k_; ++i){
		res *= (n - i)/(double)(i + 1u);
	}
	return res;
}

//Returns the number of possible combinations of swaps in a sublattice
//given species and site counts in the sublattice (returns double since I'll 
//be dividing w/ it)
double CountCombos(const ushort nSpeciesTot, const ushort* speciesArr, 
				   const ushort nSites){
	double res = 1.0;
	for(ushort i = 0; i < nSpeciesTot; ++i){
		if(speciesArr[i] == 0u){
			continue;
		}
		uint nSum = 0u;
		for(ushort j = 0; j < i; ++j){
			nSum += speciesArr[j];
		}
		res *= NCK(nSites - nSum, speciesArr[i]);
	}
	return res;
}

//Givs a uniform mapping of probabilities for random weighted choices 
//of discrete distributions
//pLen is the length of absolute probabilities p (sum p = 1.0)
//places the distribution into two arrays: lower and upper bounds
//meant to be used by choosing the only index i that satisfies:
//lower[i] <= rand([0, 1)) < upper[i] 
void UniformMap(double** lo, double** hi, const ushort pLen, const double* p){
	*lo = malloc(pLen*sizeof(double));
	*hi = malloc(pLen*sizeof(double));

	(*lo)[0] = 0.0; 
	(*hi)[0] = p[0];
	for(ushort i = 1; i < pLen; ++i){
		(*lo)[i] = (*hi)[i - 1];
		(*hi)[i] = (*hi)[i - 1] + p[i];
	}
	return;
}

//Returns a random integer chosen from a discrete distribution whose
//bounds are given and described by UniformMap()
ushort RInt(const ushort len, 
			const double*restrict lo, const double*restrict hi){
	double r = (double)rand()/(double)RAND_MAX;
	for(ushort i = 0; i < len; ++i){
		if((lo[i] <= r) && (r < hi[i]))
			return i;
	}
	///We might get here on account of floating point errors, esp. if
	///some sublattices have a "0.0" prob. of swaps.
	///In that case, 0 is always a good choice since we're guarenteed to 
	///always have at least one sublattice, and the 0th is most likley the
	///one with the highest prob. anyways
	return (ushort)0;
}
