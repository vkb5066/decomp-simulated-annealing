//Insertion sorting algorithm for doubles
//Written because often, the distance array that needs sorted is relativly 
//small in size.  You may want to use this instead of the std's qsort, which
//has a lot of overhead

#include <stdlib.h>
#include <string.h>
#include "Constants.h"


void ISort_d(double** arr, const ushort len){
	double* a; double* b; double* tmp = malloc(sizeof(double));
	for(ushort i = 0; i < len; ++i){
		for(ushort j = i; j > (ushort)0;){
			b = *arr + j;
			a = *arr + --j;
			if(*b < *a){
				memcpy(tmp, a, sizeof(double));
				memcpy(a, b, sizeof(double));
				memcpy(b, tmp, sizeof(double));
			}
		}
	}
	free(tmp);
	return;
}