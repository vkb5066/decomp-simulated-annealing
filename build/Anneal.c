#include <stdlib.h>
#include <math.h>

#include "Decs.h"
#include "Constants.h"

//Controls swapping, annealing, etc

//Tries to make a random, meaningful swap between two species in a sublattice
//Stores the sucessful swap indices in spInd1, spInd2
void RSwap(sublat*restrict*restrict s, 
		   ushort*restrict spInd1, ushort*restrict spInd2){
	for(uint i = 0; i < N_GIVEUP; ++i){
		int n1 = rand() % ((*s)->nSites);
		int n2 = rand() % ((*s)->nSites);
		if(*(*s)->sites[n1]->species != *(*s)->sites[n2]->species){
			ushort tmp = *(*s)->sites[n1]->species;
			*(*s)->sites[n1]->species = *(*s)->sites[n2]->species;
			*(*s)->sites[n2]->species = tmp;
			*spInd1 = (ushort)n1;
			*spInd2 = (ushort)n2;
			return;
		}
	}
	///We get here on some rare occasions.  On even rarer occasions, not
	//updating spInd* will result in crashes (the last sucessful swap was on
	//a sublattice with 64 sites so spInd1 = 63, spInd2 = 42.  Now, this sublat
	//only has 32 sites, but we've exceeded N_GIVEUP so spInd* aren't updated.
	//So when we try to index the new 32-site sublattice with the old spInds, 
	//there is an OOB error.
	//Each sublattice should always have at least one site, so 0 is a good 
	//fallback option
	*spInd1 = (ushort)0;
	*spInd2 = (ushort)0;
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
void UniformRandSwaps(const ushort nSublats, 
					  sublat*restrict*restrict*restrict sublats,
					  const double*restrict loBounds, 
					  const double*restrict hiBounds, 
					  const uint nSwaps,
					  ushort*restrict slInd, 
					  ushort*restrict spInd1, ushort*restrict spInd2){
	for(uint i = 0; i < nSwaps; ++i){
		*slInd = RInt(nSublats, loBounds, hiBounds);
		RSwap(&(*sublats)[*slInd], &(*spInd1), &(*spInd2));
	}
	return;
}

//Takes the just recently swapped sublattice indices, and the index of the
//sites in that index, and marks the corresponding environments as needing 
//their comparison arrays recalculated
void MarkEnvsToRecalc(const sublat** sublats,
					  const ushort swapSublatInd, const ushort swapSiteInd1,
					  const ushort swapSiteInd2){
	for(ushort i = 0; 
		i < sublats[swapSublatInd]->sites[swapSiteInd1]->nParentEnvs; ++i){
		sublats[swapSublatInd]->sites[swapSiteInd1]->parentEnvs[i]->recalc = 
																	 (ushort)1;
	}
	for(ushort i = 0;
		i < sublats[swapSublatInd]->sites[swapSiteInd2]->nParentEnvs; ++i){
		sublats[swapSublatInd]->sites[swapSiteInd2]->parentEnvs[i]->recalc =
																	 (ushort)1;
	}
}


//Predicts the energy of a configuration by matching it's environments (trial) 
//to environments of already known energies (known).  Returns the total energy
//and, in the case that there are unknown energies, populates an array 
//(initially all 0s) with a 1 for every unknown entry i.e. {0 0 1 0 0 1} means
//that trialEnv[2] and trialEnv[5] have no matches in the known array
void UpdateNrg(double*restrict lastNrg, const ushort nSites,
			   const hbkt*restrict table,
			   const ushort nTrialEnvs, 
			   env*restrict*restrict*restrict trialEnvs,
			   ushort*restrict nUnknowns, ushort*restrict*restrict unknowns,
			   const ushort nSpecTot,
			   ushort*restrict nRecalcInds, 
			   ushort*restrict*restrict recalcInds){
	(*nUnknowns) = (ushort)0;
	(*nRecalcInds) = (ushort)0;

	ushort found;
	//For each trial env (that needs recalculated)...
	(*lastNrg) *= (double)nSites;
	for(ushort i = 0; i < nTrialEnvs; ++i){
		if((*trialEnvs)[i]->recalc){
			SetSpecArr_e(&((*trialEnvs)[i]), nSpecTot); 
			SetDistArr_e(&((*trialEnvs)[i]));
		
			///Unmark for recalculation
			(*trialEnvs)[i]->recalc = (ushort)0;
			///Mark as having been recalculated
			(*recalcInds)[*nRecalcInds] = i;
			(*nRecalcInds)++;

			///...find corresponding fit env and ajdust energy accordingly
			env* fit = SearchEnv(table, (const env*)(*trialEnvs)[i], nSpecTot);
			if(fit){
				(*lastNrg) -= (*trialEnvs)[i]->nrg; ///subtract prev nrg
				(*trialEnvs)[i]->nrg = fit->nrg; ///set new nrg
				(*lastNrg) += fit->nrg; ///add new nrg
			}
			else{
				(*nUnknowns)++;
				(*unknowns)[i] = 1;
#if ASSIGN_UNK_ENV_NRG
				(*lastNrg) += UNK_ENV_NRG;
#endif
			}
		}
	}
	(*lastNrg) /= (double)nSites;

	return;
}

//step: current step, out of nStepsTot
//tau*: thermal energy beginning and ending point, eV / site
double CalcTau(const uvlong step, const uvlong nStepsTot, 
			   const double tauStart){
	double fac = 1.0 - ((double)step/(double)nStepsTot);
	return tauStart*fac*fac*fac*fac;
}

//Deals with organizing, adding to, or overwriting the array of best
//sites and their energies
//increments a success counter for each sucessful addition
void BestSiteHandler(double*restrict*restrict bstNrgs, 
					 ushort*restrict*restrict*restrict bstOccs, 
					 const uint bstLen, 
					 const site*restrict*restrict allSites,
					 const ushort nSites,
					 const double currNrg, uint*restrict suc){
	//Early termination, hit most often:
	//If the current energy is higher than any of the 
	//energies in the array, don't bother doing the expensive searches
	if((*bstNrgs)[bstLen - 1u] < currNrg) return;
	
	//Check to make sure this structure is unique: if not, skip
	//Note: this might be better as a hash search...but theoretically we
	//don't get here very often, so it might be fine
	for(uint i = 0; i < bstLen; ++i){
		for(ushort j = 0; j < nSites; ++j){
			if((*bstOccs)[i][j] != *(allSites[j]->species)){
				goto Next;
			}
		}
		return; ///if we hit this, we've found a duplicate entry
		Next: NOP
	}
	
	//Get index to replace
	uint repInd = SearchB_d(*bstNrgs, bstLen, currNrg);

	//Shift the best energies and best occupancies to make room for new
	//one if necessary
	Insert_d(&(*bstNrgs), bstLen, repInd, currNrg);
	Insert_o(&(*bstOccs), bstLen, repInd, allSites, nSites);

	(*suc)++;
	return;
}
