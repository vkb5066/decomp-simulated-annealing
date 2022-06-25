//General
typedef unsigned int uint;
typedef unsigned short ushort;
typedef struct site site;
typedef struct env env;
typedef struct sublat sublat;
typedef struct decomp decomp;

#define NOP ;

//Input-Output
#define BASE 10u
#define INFILE_LINE_MAX 1024u ///lines to read in infile before aborting
#define LINESIZE 255u ///max len of any given line in input file
#define EOL '\n'

//Size of max element ID string
#define ELEMSIZE 8u ///max number of characters in an element string

//Small element for approximate comparisons
#define DRCT_EPS 1.0E-3 ///fractional
#define CART_EPS 5.0E-4 ///cartsian

//Number of times we try to make a "significant" random swap before giving up
//By "significant", I mean that two species that are swapped are not identical
#define N_GIVEUP 128u

//Number of swaps from one structure to the next for the next to be considered
//a neighbor of the first.  Should just leave this as one, probably
#define N_NBR_SWAPS 1u

//Number of unique enviornments added to the (large) array before 
//we attempt to  re-allocate memory with an additional N_REALLOC
//spaces
#define N_REALLOC 256u