#include <stdlib.h>

#include "Decs.h"
#include "Constants.h"

//Controls swapping, annealing, etc

//Tries to make a random, meaningful swap between two species in a sublattice
//Stores the sucessful swap indices in spInd1, spInd2
void RSwap(sublat** s, ushort* spInd1, ushort* spInd2){
	for(uint i = 0; i < N_GIVEUP; ++i){
		int n1 = rand() % ((*s)->nSites);
		int n2 = rand() % ((*s)->nSites);
		if(*(*s)->sites[n1]->species != *(*s)->sites[n2]->species){
			ushort tmp = *(*s)->sites[n1]->species;
			*(*s)->sites[n1]->species = *(*s)->sites[n2]->species;
			*(*s)->sites[n2]->species = tmp;
			(*spInd1) = (ushort)n1;
			(*spInd2) = (ushort)n2;
			return;
		}
	}
	return;
}

//For first initialization only.  Random swaps over all sublattices
//disregarding probability distributions (unless it is super low, in which
//case we don't even try)
void NonUniformRandSwaps(const ushort nSublats, sublat*** sublats){
	ushort dummy1 = 0; ushort dummy2 = 0;
	for(ushort i = 0; i < nSublats; ++i){
		for(ushort j = 0; j < (*sublats)[i]->nSites; ++j) {
			RSwap(&(*sublats)[i], &dummy1, &dummy2);
		}
	}
	return;
}

//Attempts to preform nSwaps random (important) swaps by choosing
//(first) a sublattice to consider - weighted distribution
//(second) two sites in the sublattice s.t. the species are unequal
//         assuming a uniform distribution of swap probabilities in this case
//Adds the sublattice index, index 1 and index 2 of the last swapped elements
void UniformRandSwaps(const ushort nSublats, sublat*** sublats,
					  const double* loBounds, const double* hiBounds, 
					  const uint nSwaps,
					  ushort* slInd, ushort* spInd1, ushort* spInd2){
	for(uint i = 0; i < nSwaps; ++i){
		*slInd = RInt(nSublats, loBounds, hiBounds);
		RSwap(&(*sublats)[*slInd], &(*spInd1), &(*spInd2));
	}
	return;
}

//Predicts the energy of a configuration by matching it's environments (trial) 
//to environments of already known energies (known).  Returns the total energy
//and, in the case that there are unknown energies, populates an array with the
//(initially all 0s) with a 1 for every unknown entry i.e. {0 0 1 0 0 1} means
//that trialEnv[2] and trialEnv[5] have no matches in the known array
double CalcNrg(const uint nKnownEnvs, const env** knownEnvs, 
			   const ushort nTrialEnvs, env*** trialEnvs,
			   ushort* nUnknowns, ushort** unknowns,
			   const ushort nSpecTot){
	(*nUnknowns) = 0;

	double nrg = 0.0;
	ushort found;
	for(ushort i = 0; i < nTrialEnvs; ++i){
		SetSpecArr_e(&((*trialEnvs)[i]), nSpecTot); 
		SetDistArr_e(&((*trialEnvs)[i]));
		
		found = (ushort)0;
		for(uint j = 0; j < nKnownEnvs; ++j){
			if(Cmpr_e((*trialEnvs)[i], knownEnvs[j], nSpecTot)){
				nrg += knownEnvs[j]->nrg;
				found = 1u;
				break;
				
			}
		}

		if(!found){
			(*nUnknowns)++;
			(*unknowns)[i] = 1;
		}
	}

	return nrg;
}

//step: current step, out of nStepsTot
//tau*: thermal energy beginning and ending point, eV / site
//cool*: fraction (0 -> 1) of time before cooling begins, ends
//       coolBegin < coolEnd
//alpha: controls the shape of the cooling schedule curve.
//       reasonable values are between 2/3 and 5 (higher = more rapid cooling)
double CalcTau(const uint step, const uint nStepsTot, 
			   const double tauStart, const double tauEnd,
			   const double coolBegin, const double coolEnd,
			   const double alpha){
	double scx = ((double)step/(double)nStepsTot - coolBegin) / 
				 (coolEnd - coolBegin); 
	return (tauStart + (tauEnd - tauStart)*SwitchFunction(alpha, scx));
}
