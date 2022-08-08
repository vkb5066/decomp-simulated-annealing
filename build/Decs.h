#include "Constants.h"


//Plain Functions--------------------------------------------------------------
//File Input/output (IO.c)
void ReadLatticeFile(const char* fileName, 
					 double A[3][3], 
					 ushort* nSpecies, 
					 ushort* nSites, site*** sites,
					 ushort* nSublats, sublat*** sublats);
ushort ReadParamFile(const char* fileName, double* rCut,
					 uvlong* nSteps, uvlong* startStep,
					 double* startNrg, int* rSeed,
					 ushort** initState);
void ReadEnvFile(const char* fileName, 
				 uint* nEnvs, env*** envs,
				 const ushort nSpecTot, const ushort* fixSpecArr);
void ReadWarmStartFile(const char* fileName,
					   uint* cSizeE, uint* mSizeE, env*** arrE,
					   uint* cSizeD, uint* mSizeD, uint*** arrD,
					   ushort*** arrS,
					   uint* nReppedEnvs,
					   const ushort nSpecTot, const ushort* fixSpecArr,
					   const ushort nSites);

void UnknownEnvHandler(const ushort nSites, const site** sites,
					   const uvlong step,
					   const ushort nUnknowns, const ushort* unknownEnvs,
					   const ushort nTrialEnvs, const env** trialEnvs);
void SitePrintout(const uint nPrints, const ushort nSites, const ushort** arr);
void EnergyPrintout(const uint nPrints, const double* arr);
void DecPrintout(const uint nPrints, const uint maxCol,
				 const ushort nSites, const uint** arr);
void EnvPrintout(const uint nEnvs, const env** envs,
				 ushort giveNrg, ushort giveRep);

//Lattice functions
void MoveToCell(site** s);
void GiveNImgs(ushort nImgs[3], const double rCut, const double A[3][3]);
void SetSiteGeoms(const ushort nEnvs, env*** envs,
				  const ushort nSites, site*** sites,
				  const ushort nImgs[3], const double rCut,
				  const double A[3][3], const ushort nSpecTot);
void SetConfig(const ushort nSites, site*** sites, const ushort* species);
void DecompConv(const ushort nTotEnvs, const ushort* full, uint** sparse);
int DecompCmpr_fs(const ushort nTotEnvs,
				  const ushort* full, const uint* sparse);
int DecompCmpr_ss(const ushort nTotEnvs, const uint* a, const uint* b);
void DeepcopyState(const ushort nSites, const site** src, site*** dst);

//Math
double NCK(const uint n, const uint k);
double CountCombos(const ushort nSpeciesTot, const ushort* speciesArr,
				   const ushort nSites);
void UniformMap(double** lo, double** hi, const ushort pLen, const double* p);
ushort RInt(const ushort len, 
			const double*restrict lo, const double*restrict hi);

//Annealing
void RSwap(sublat*restrict*restrict s, 
		   ushort*restrict spInd1, ushort*restrict spInd2);
void NonUniformRandSwaps(const ushort nSublats, sublat*** sublats);
void UniformRandSwaps(const ushort nSublats, 
					  sublat*restrict*restrict*restrict sublats,
					  const double*restrict loBounds, 
					  const double*restrict hiBounds,
					  const uint nSwaps,
					  ushort*restrict slInd, 
					  ushort*restrict spInd1, ushort*restrict spInd2);
void MarkEnvsToRecalc(const sublat** sublats,
					  const ushort swapSublatInd, const ushort swapSiteInd1,
					  const ushort swapSiteInd2);
void UpdateNrg(double* restrict lastNrg, const ushort nSites,
			   const hbkt* restrict table,
			   const ushort nTrialEnvs,
			   env* restrict* restrict* restrict trialEnvs,
			   ushort* restrict nUnknowns, ushort* restrict* restrict unknowns,
			   const ushort nSpecTot,
			   ushort* restrict nRecalcInds,
			   ushort* restrict* restrict recalcInds);
double CalcTau(const uvlong step, const uvlong nStepsTot,
			   const double tauStart);
void BestSiteHandler(double*restrict*restrict bstNrgs,
					 ushort*restrict*restrict*restrict bstOccs,
					 const uint bstLen,
					 const site*restrict*restrict allSites, 
					 const ushort nSites,
					 const double currNrg, uint*restrict suc);

//Sorting, searching
void ISort_d(double** arr, const ushort len);
uint SearchB_d(const double* arr, const uint len, const double val);
void Insert_d(double** arr, const uint len, const uint ind, const double val);
void Insert_o(ushort*restrict*restrict*restrict arr,
			  const uint len, const uint ind,
			  const site*restrict*restrict val, const ushort valLen);


//-----------------------------------------------------------------------------


//Structs----------------------------------------------------------------------
//Site structure
struct site{
	ushort* species; ///identifier for this element
	double crdsD[3]; ///direct coordinates a, b, c
	double crdsC[3]; ///cartesian coordinates x, y, z

	ushort nParentEnvs;
	env** parentEnvs;

	site* self; ///points to this site's original site
};
void SetCarCrds_s(site** s, const double A[3][3]);
double DSqd_s(const site* a, const site* b);
void Print_s(const site* s);

//Environment structure (holds sites)
struct env{
	ushort nSites;
	site** sites;

	ushort*restrict sortedSpecies;
	ushort nDists;
	double*restrict sortedDists;

	double nrg;

	uint decId;
	ushort repped;
	ushort recalc;
};
void AllocCmpArrs_e(env** e, const ushort nSpecTot);
void SetSpecArr_e(env** e, const ushort nSpecTot);
void SetDistArr_e(env** e);
int Cmpr_e(const env* e1, const env* e2, const ushort nSpecTot);
void Print_e(const env* e, ushort giveNrg, ushort giveRep);
void SwapCmpArrs_e(const env*restrict src, env*restrict*restrict dst, 
				   const ushort nSpecTot);

//Sublattice structure (holds sites)
struct sublat{
	ushort nSites;
	site** sites;

	double swapProb;
};

//Search tree structure (holds pointers to envs)
struct branch{
	ushort nBranches;
	branch** branches;

	double key;
	env* thisEnv;
};
void Alloc_b(branch** b, double key);
void Dealloc_b(branch** b);
env* Add_b(branch*restrict*restrict head, const env*restrict toAdd, 
		   uint*restrict success);
env* Search_b(const branch*restrict head, const env*restrict searchKey);

//Hash table functions (holds pointers to trees)
struct hbkt{
	ushort nEntries;
	ushort** specArrs;
	branch** branchHeads;
};
ushort ArrToHash_h(const ushort len, const ushort* arr);
void InitTable_h(hbkt* ptr);
void DeallocTable_h(hbkt* h);
branch* Add_h(hbkt*restrict table, const ushort hashInd,
			  const ushort specArrLen, const ushort*restrict specArr,
			  ushort*restrict collisions);
branch* Search_h(const hbkt*restrict table, const ushort hashInd,
				 const ushort specArrLen, const ushort*restrict specArr);

struct gbkt{
	uint nEntries;
	uint** sparseArrs;
};
ushort ArrToHash_gh(const ushort len, const ushort* arr);
void InitTable_gh(gbkt* ptr);
void DeallocTable_gh(gbkt* g);
void Add_gh(gbkt*restrict table, const uint*restrict sparseArr,
			const ushort nTotEnvs, uint*restrict colls,
			ushort*restrict hash, uint*restrict addInd, ushort*restrict suc);

struct obkt{
	uint nEntries;
	ushort** occArrs;
};
ushort ArrToHash_oh(const ushort len, const ushort* arr);
void InitTable_oh(obkt* ptr);
void DeallocTable_oh(obkt* o);
void Add_oh(obkt*restrict table, const ushort*restrict sites,
			const ushort nSites, ushort*restrict suc);

//General add, search functions
env* AddEnv(hbkt*restrict table, const ushort nSpecTot, 
			const env*restrict toAdd,
			ushort*restrict collCnt, uint*restrict newEnvCount);
env* SearchEnv(const hbkt*restrict table, const env*restrict query, 
			   const ushort nSpecTot);

//-----------------------------------------------------------------------------