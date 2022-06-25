#include "Constants.h"


//Plain Functions--------------------------------------------------------------
//File Input/output (IO.c)
void ReadLatticeFile(const char* fileName, 
					 double A[3][3], 
					 ushort* nSpecies, 
					 ushort* nSites, site*** sites,
					 ushort* nSublats, sublat*** sublats);
ushort ReadParamFile(const char* fileName, double* rCut,
					 uint* nSteps, uint* startStep,
					 double* startNrg, double* endNrg,
					 double* startCool, double* endCool,
					 double* alpha, int* rSeed,
					 ushort** initState);
void ReadEnvFile(const char* fileName, 
				 uint* nEnvs, env*** envs,
				 const ushort nSpecTot, const ushort* fixSpecArr);
void ReadWarmStartFile(const char* fileName,
					   uint* cSizeE, uint* mSizeE, env*** arrE,
					   uint* cSizeD, uint* mSizeD, uint*** arrD,
					   ushort*** arrS,
					   uint* unreppedEnvCount,
					   const ushort nSpecTot, const ushort* fixSpecArr,
					   const ushort nSites);

void UnknownEnvHandler(const ushort nSites, const site** sites,
					   const uint step,
					   const ushort nUnknowns, const ushort* unknownEnvs,
					   const ushort nTrialEnvs, const env** trialEnvs);
void SitePrintout(const uint nPrints, const ushort nSites, const ushort** arr);
void DecPrintout(const uint nPrints, const uint maxCol,
				 const ushort nSites, const uint** arr);
void EnvPrintout(const uint nEnvs, const env** envs,
				 ushort giveNrg, ushort giveRep);

//Lattice functions
void MoveToCell(site** s);
void GiveNImgs(ushort nImgs[3], const double rCut, const double A[3][3]);
void SetSiteGeoms(const ushort nEnvs, env*** envs,
				  const ushort nSites, const site** sites,
				  const ushort nImgs[3], const double rCut,
				  const double A[3][3], const ushort nSpecTot);
void SetConfig(const ushort nSites, site*** sites, const ushort* species);
void DecompConv(const ushort nTotEnvs, const ushort* full, uint** sparse);
int DecompCmpr_fs(const ushort nTotEnvs,
				  const ushort* full, const uint* sparse);

//Math
double NCK(const uint n, const uint k);
double CountCombos(const ushort nSpeciesTot, const ushort* speciesArr,
				   const ushort nSites);
void UniformMap(double** lo, double** hi, const ushort pLen, const double* p);
ushort RInt(const ushort len, const double* lo, const double* hi);
double Smooth(const double alpha, const double x);
double SwitchFunction(const double alpha, const double x);

//Annealing
void RSwap(sublat** s, ushort* spInd1, ushort* spInd2);
void NonUniformRandSwaps(const ushort nSublats, sublat*** sublats);
void UniformRandSwaps(const ushort nSublats, sublat*** sublats,
					  const double* loBounds, const double* hiBounds,
					  const uint nSwaps,
					  ushort* slInd, ushort* spInd1, ushort* spInd2);
double CalcNrg(const uint nKnownEnvs, const env** knownEnvs,
			   const ushort nTrialEnvs, env*** trialEnvs,
			   ushort* nUnknowns, ushort** unknowns,
			   const ushort nSpecTot);
double CalcTau(const uint step, const uint nStepsTot,
			   const double tauStart, const double tauEnd,
			   const double coolBegin, const double coolEnd,
			   const double alpha);

//-----------------------------------------------------------------------------


//Structs----------------------------------------------------------------------
//Site structure
struct site{
	ushort* species; ///identifier for this element
	double crdsD[3]; ///direct coordinates a, b, c
	double crdsC[3]; ///cartesian coordinates x, y, z

	site* self; ///points to this site's original site
};
void SetCarCrds_s(site** s, const double A[3][3]);
double DSqd_s(const site* a, const site* b);
void Print_s(const site* s);

//Environment structure (holds sites)
struct env{
	ushort nSites;
	site** sites;

	ushort* sortedSpecies;
	ushort nDists;
	double* sortedDists;

	double nrg;

	ushort repped;
};
void AllocCmpArrs_e(env** e, const ushort nSpecTot);
void SetSpecArr_e(env** e, const ushort nSpecTot);
void SetDistArr_e(env** e);
int Cmpr_e(const env* e1, const env* e2, const ushort nSpecTot);
void Print_e(const env* e, ushort giveNrg, ushort giveRep);

//Sublattice structure (holds sites)
struct sublat{
	ushort nSites;
	site** sites;

	double swapProb;
};

//-----------------------------------------------------------------------------