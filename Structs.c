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
	qsort((*e)->sortedDists, (*e)->nDists, sizeof(double), dblcmpr);
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