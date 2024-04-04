/****************************************************************************/
/* These are the library functions that are used within the calculation of the band structure */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <mkl.h>
#include <sys/time.h>
#include <time.h>
#include <assert.h>
#include "complex.h"
#include <omp.h>

/****************************************************************************/
/* These are some common unit conversion and multiplication schemes that help make the program more readable */

#define sqr(x)       ((x) * (x))
#define cube(x)      ((x) * (x) * (x))
#define forth(x)     ((x) * (x) * (x) * (x))
#define ISWAP(a,b)   {double temp = (a); (a) = (b); (b) = -temp;}
#define AUTOEV    27.2114
#define AUTONM	  0.05291772108
#define PIE       3.14159265358979323846
#define TWOPI     6.28318530717958647692
#define HBAR	  1.0
#define NDIM 	  3
#define MASS 	  1.0
#define EPS     1e-10
#define NPParams 6
//#define LOG     "fit.log"

/****************************************************************************/
/* Structure declarations */

/*typedef struct complexnumber {
  double re, im;
} complexnumber;
*/
typedef struct vector {
  double x, y, z;
  double mag;
} vector;

typedef struct atom {
  char symbol[1000];
  int type;
  vector pos;
  vector defPos;
  double ppParams[NPParams];
  double spinOrbit;
  double fF;
} atom;

typedef struct param_ {
  int nBandStructures, nIterations;
  int nAtoms, nAtomTypes, nBands;
  int nKPoints, nKPE; 
  int nkp[20], ne[20], na[20];
  int nqp[20], nbv[20], nmbv[20], nElec[20];
  double s[20], mke[20], vol[20];
  double defS[20], defVol[20]; // DJ
  int nQpoints, nBasisVectors, nMaxBV;
  double scale, maxKE, beta, stepSize, stepSizeScale[NPParams], defWeight[4];
  int nThreads;
  double SOWidth;
} param;

/****************************************************************************/
/* Function declarations */

/* Functions that initialize the system */
/* read.c */
void readInput(param *params);
void readMultiInput(param *params);
void copyBaseParamsToArrays(param *params);
void readSystem(param *params, vector *unitCellvectors, vector *defUnitCellvectors, atom *atoms, int fileIndex);
void readExpBandStructure(double *expBandStructure, double *expDefPot, int nKPoints, int nElectrons, int fileIndex);
void readKPoints(vector *kPoint, vector *gVectors,int nKPoints, int fileIndex);

/* init.c */
int newAtomType(atom *atoms, int i);
double getCellVolume(vector *cellVectors);
void getRecipSpaceBasisVectors(vector *unitCellvectors, vector *gVectors, param *params, int fileIndex);
void initAtomicFormFactors(param params, atom *atoms);

/* states.c */
void calcBasisVectors(vector *tmpBV, vector *gVectors, param *params, int fileIndex);
void finalizeBasisStates(vector *basisStates, vector *tmpBV, param *params, int fileIndex);
int countDistinctMagnitudes(vector *vect, int nVectors);
void quickSortVectors(vector *vectors, int l, int r);
int partitionVectors(vector *vectors, int l, int r);

/* Functions that involve the Hamiltonian matrix */
/* hamiltonian.c */
void calcHamiltonianMatrix(complexnumber *hamiltonian, complexnumber *localPotential, complexnumber* spinOrbitMatrix, vector *basisStates, vector *kVector, int n, int SOBool);
void calcPotentialEnergyMatrix(complexnumber *localPotential, vector *basisStates, atom *atoms, int nbv, int nAtoms, double cellVolume, double defCellVolume, int defBool);
double calcPot(double q, double *param);
double calcStrainPot(double q, double *param, double vol, double defVol);
void zeroComplexMatrix(complexnumber *mat, int n1, int n2);
void calcVSpinOrbit(double *vSpinOrbit, vector *basisStates, vector *kPoint, param params, int ibs); 
void calcSpinOrbitPotentialMatrix(complexnumber *spinOrbitMatrix, double *vSpinOrbit, vector *basisStates, vector kVector, int kindex, atom *atoms, param params,int ibs);

void calcKineticEnergyMatrix(complexnumber *hamiltonian, vector *basisStates, vector kVector, param params);

/* diag.c */
void diagonalizeHamiltonian(complexnumber *hamiltonian, double *energies, int nBV);


double calcBessel(double *bessely, double *besselm, double splwdth, double splmax, double x);
void calcCohenPEM(complexnumber *localPotential, vector *basisStates, atom *atoms, param params, vector* unitCellvectors);
complexnumber cohenLookUp(double gSqrMag, param params);
void calcZungerVsVaPEM(complexnumber *localPotential, vector *basisStates, atom *atoms, param params);
complexnumber zungerLookUp(double gMag, param params);
void calcRealSpacePotential(atom *atoms, param params);

/* Functions needed for the Monte Carlo fitting process */
/* fit.c */
void runMonteCarlo(complexnumber *hamiltonian, complexnumber *potential, complexnumber *spinOrbitPot, double *vSpinOrbit, vector *basisStates, 
	vector *defBasisStates, double *expBandStructure, double *expDefPot, vector *kPoint, double *bandStructure, double *energies, atom *atoms,
	param params, int SOBool, vector *gVectors);
void calcBandStructure(double *bandStructure, complexnumber *hamiltonian, complexnumber *potential, complexnumber *spinOrbitPot, double *vSpinOrbit,
        vector *basisStates, vector *kPoint, double *energies, atom *atoms, param params, int SOBool);
void calcBandGaps(double *cBands, double *vBands, complexnumber *hamiltonian, complexnumber *potential, complexnumber *spinOrbitPot,
        vector *basisStates, double *energies, atom *atoms, param params, double *cellVolumes, double *defCellVolumes, int defBool);
void calcDeformationPotentials(double *defPot, complexnumber *hamiltonian, complexnumber *potential, complexnumber *spinOrbtiPot,
        vector *basisStates, vector *defBasisStates, double *energies, atom *atoms, param params);
double calcAllIndividualWeightedMSE(double *bandStructure, double *expBandStructure, param params);
double calcWeightedMSE(double *array1, double *array2, double *weights, int nElectrons, int nKPoints);
double calcAllIndividualMSE(double *bandStructure, double *expBandStructure, vector *kPoint, double *defPot, double *expDefPot, atom *atoms, param params, vector *gVectors);
double calcMSE(double *array1, double *array2, vector *kVector, int arrayLength, vector *gVectors);
void savePseudoParams(double *tmpPseudos, atom *atoms, param params);
void updatePseudoParams(atom *atoms, param params);
void copyPreviousAtomParams(atom *atoms, int i);
void revertPseudoParams(atom *atoms, double *tmpPseudos, param params);
void storeEnergies(double *bandStructure, double *energies, int index, int length);
double randNumber(void);

/* Functions that write the output */
void writeResults(double *bandStructure, vector *kPoint, double bestMSE, atom *atoms, param params);
void writeIteration(double MSE, int iteration, atom *atoms, param params, double *defPot, FILE *pf);
void writeInput(atom *atoms, param params, FILE *pf);
void writePseudoParams(atom *atoms, param params, FILE *pf);
void writeBestPseudoParams(atom *atoms, param params);
void writeVector(vector vect, FILE *pf);
void writeCurrentTime(FILE *pf);
void printCurrentTime();
void writeSeparation(FILE *pf);

/* Functions for vector structures */
vector retAddedVectors(vector vect1, vector vect2);
vector retSubtractedVectors(vector vect1, vector vect2);
vector retScaledVector(vector vect, double scale);
vector retCrossProduct(vector vect1, vector vect2);
vector retZeroVector(void);
double retDotProduct(vector vect1, vector vect2);
int compareVectorMagnitudes(vector vect1, vector vect2);

/****************************************************************************/
