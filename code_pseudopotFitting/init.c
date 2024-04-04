/****************************************************************************/
/* This files does the required initialization before calculating the hamiltonian */

#include "fit.h"

/****************************************************************************/

int newAtomType(atom *atoms, int i) {
  int j;
 
  for (j = 0; j < i; j++) {
    if (! strcmp(atoms[j].symbol, atoms[i].symbol)) return 0;
  }
 
  return 1;
}

/****************************************************************************/

double getCellVolume(vector *cellVectors) {
  double cellVolume;
  vector temp;
  
  temp = retCrossProduct(cellVectors[1], cellVectors[2]);
  cellVolume = retDotProduct(cellVectors[0], temp);
  
  return (cellVolume);
}

/****************************************************************************/

void getRecipSpaceBasisVectors(vector *unitCellvectors, vector *gVectors, param *params, int fileIndex) {
  double preFactor = TWOPI / params->vol[fileIndex];
  fflush(stdout);

  gVectors[0] = retCrossProduct(unitCellvectors[1], unitCellvectors[2]);
  gVectors[0] = retScaledVector(gVectors[0], preFactor);
  gVectors[1] = retCrossProduct(unitCellvectors[2], unitCellvectors[0]);
  gVectors[1] = retScaledVector(gVectors[1], preFactor);
  gVectors[2] = retCrossProduct(unitCellvectors[0], unitCellvectors[1]);
  gVectors[2] = retScaledVector(gVectors[2], preFactor); 
  

  //TODO: why does this happen here? move this somwhere more sane?
  /* overestimate of the number of basis vectors */
  double minGMag = 100000.0; int i;
  
  for (i = 0; i < NDIM; i++) {
      if (gVectors[i].mag < minGMag) minGMag = gVectors[i].mag;
  }
    
  params->nmbv[fileIndex] = (int) (sqrt(2*params->mke[fileIndex]) / minGMag);

  //printf("BVMaxN=%d\n",params->BVMaxN);
  return;
}

/****************************************************************************/

void initAtomicFormFactors(param params, atom *atoms) {
  FILE *pf;
  int i, j, k;
  char fileName[1000];

  int ina = 0;
  for (k = 0; k < params.nBandStructures; k++) {
    for (i = 0; i < params.na[k]; i++) {
      if (newAtomType(atoms, ina)) {
        memset(&fileName[0], 0, sizeof(fileName));
        strcpy(fileName, atoms[ina].symbol);
        strcat(fileName, "Params.par");
        pf = fopen(fileName, "r");
        for (j = 0; j < NPParams+1; j++) {
          if (j<NPParams) {
              fscanf(pf, "%lg", &atoms[ina].ppParams[j]);
          }
          else {
	    fscanf(pf, "%lg", &atoms[ina].spinOrbit);
	  }
        }
        fclose(pf);
      }
      else {
        copyPreviousAtomParams(atoms, ina);
      }
      ina++;
    }
  }
  writeSeparation(stdout);
  fprintf(stdout, "Initial pseudopotentials Params:\n");
  writePseudoParams(atoms, params, stdout);
  fprintf(stdout, "\n" );
  fflush(stdout);
  return;
}

