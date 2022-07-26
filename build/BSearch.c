
#include <stdlib.h>
#include <string.h>

#include "Decs.h"
#include "Constants.h"


//Searches a sorted array and returns the index i of the array that satisfies
//arr[i - 1] <= val < arr[i + 1]
//i.e. the returned value is the index to insert val into to ensure that the
//array is sorted after an insertion
uint SearchB_d(const double* arr, const uint len, const double val){
	if(val > arr[len - 1u]) return len; ///most likley case
	if(val < arr[0]) return 0u; ///least likley case, but keep this here

	uint lo = 0u; uint hi = len - 1u; uint mi;
	while(lo <= hi){
		mi = (lo + hi)/2u;

		if((arr[mi] <= val) && (val <= arr[mi + 1u])){
			return mi + 1u;
		}
			
		if(arr[mi] < val)  lo = mi + 1u;
		else hi = mi - 1u;
	}

	return (uint)-1; ///ruh roh
}

//Inserts a double into a sorted array to maintain it's order by shifting
//all elements @ + above ind up by one
//OVERWRITES the value at arr[len - 1]
void Insert_d(double** arr, const uint len, const uint ind, const double val){
	///shift
	for(uint i = len - 1u; i > ind; --i){
		//memcpy(*arr + i, *arr + i - 1u, sizeof(double));
		(*arr)[i] = (*arr)[i - 1u];
	}
	///insert
	(*arr)[ind] = val;

	return;
}

//Ditto above, but does so with the occupancy array
//OVERWRITES the value at arr[len - 1]
void Insert_o(ushort*restrict*restrict*restrict arr, 
			  const uint len, const uint ind, 
			  const site*restrict*restrict val, const ushort valLen){
	///shift
	for(uint i = len - 1u; i > ind; --i){
		for(ushort j = 0; j < valLen; ++j){
			(*arr)[i][j] = (*arr)[i - 1u][j];
		}
	}
	///insert
	for(ushort j = 0; j < valLen; ++j){
		(*arr)[ind][j] = *val[j]->species;
	}

	return;
}