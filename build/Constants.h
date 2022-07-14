//General
typedef unsigned int uint;
typedef unsigned short ushort;
typedef unsigned long long uvlong;

#ifdef _WIN32
#define restrict __restrict ///no, typdefs don't work here...
#endif

typedef struct site site;
typedef struct env env;
typedef struct sublat sublat;
typedef struct branch branch;
typedef struct hbkt hbkt;
typedef struct gbkt gbkt;

#define NOP ;

//Input-Output
#define BASE 10u
#define INFILE_LINE_MAX 1024u ///lines to read in infile before aborting
#define LINESIZE 255u ///max len of any given line in input file
#define GEN_PRINTEVERY 2048llu ///print out update every PRINTEVERY lines
#define ANN_PRINTEVERY 2048llu ///ditto ^
#define GEN_PRINTSTYLE 1 //0: no print, 1: a ton of values, 2: % done
#define ANN_PRINTSTYLE 1 ///ditto ^.  1: new lines, 2: no new lines

//Size of max element ID string
#define ELEMSIZE 8u ///max number of characters in an element string

//Small element for approximate comparisons
#define DRCT_EPS 1.0E-3 ///fractional
#define CART_EPS 5.0E-4 ///cartesian

//Whether to use std's qsort (0) or basic insertion sort (1)
//Use insertion sort for small distance tables (< 15 or so)
#define SORT_TYPE 1

//Number of times we try to make a "significant" random swap before giving up
//By "significant", I mean that two species that are swapped are not identical
#define N_GIVEUP 128u

//How to deal with unknown environments in the simulated annealing routine:
//0 (default) freaks out, won't continue the optimization, and prints out the
//offending environments
//>= 1 will assume that the unknown environment has an energy of the #define
//below it and continue the optimization
#define ASSIGN_UNK_ENV_NRG 0
#define UNK_ENV_NRG 100.0 ///eV TOTAL (not eV/site)!

//Number of unique enviornments added to arrays before attempting 
//to re-allocate memory with an additional N_REALLOC
//spaces
#define N_REALLOC 2048u

//Length of the hash table for turning arrays of species (integers) into
//a single integer.  
#define HASH_TABLE_MUL (ushort)37 ///choose a prime (31 or 37 work well)
#define HASH_TABLE_LEN (ushort)512 ///choose a power of 2

//Ditto above, but these are for the decomposition hash table
//Mind that this may be much larger than the env table
#define GHASH_TABLE_MUL (ushort)37
#define GHASH_TABLE_LEN (ushort)4096

//If false, assumes that the hash table is completly filled
//i.e. if there is one entry in a hash bin, assume that it is the 
//ONLY possible entry
//Saves at least (num total species) integer comparisons for each lookup
//Does not apply to adding a hash entry (where this assumption is dangerous)
#define SAFE_HASH_LOOKUP 1u


