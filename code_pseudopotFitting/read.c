/****************************************************************************/

#include "fit.h"

/****************************************************************************/

void readInput(param *params) {
  FILE *pf;
  int i, j;
  i = j = 0;

  char field[1000], tmp[1000];

  // Set default

  if ((pf = fopen("input.par", "r"))==NULL)
  {
    printf("Can't find input.par... exiting\n");
    exit(EXIT_FAILURE);
  }
 
  while (1) {
    fscanf(pf, "%s", field);
    if (! strcmp(field, "maxKE")) fscanf(pf, "%s %lg", tmp, &params->maxKE);
    else if (! strcmp(field, "nBandStructures")) fscanf(pf, "%s %d", tmp, &params->nBandStructures);
    else if (! strcmp(field, "nAtoms")) fscanf(pf, "%s %d", tmp, &params->nAtoms); // total for all systems
    else if (! strcmp(field, "nBands")) fscanf(pf, "%s %d", tmp, &params->nBands); // again total for all systems
    else if (! strcmp(field, "nIterations")) fscanf(pf, "%s %d", tmp, &params->nIterations); 
    else if (! strcmp(field, "nKPoints")) fscanf(pf, "%s %d", tmp, &params->nKPoints); //again total for all systems
    else if (! strcmp(field, "beta")) fscanf(pf, "%s %lg", tmp, &params->beta);
    else if (! strcmp(field, "stepSize")) fscanf(pf, "%s %lg", tmp, &params->stepSize);
    else if (! strcmp(field, "stepSizeScale")) fscanf(pf, "%s %lg %lg %lg %lg %lg %lg", tmp, &params->stepSizeScale[0], &params->stepSizeScale[1], &params->stepSizeScale[2], &params->stepSizeScale[3], &params->stepSizeScale[4], &params->stepSizeScale[5]);
    else if (! strcmp(field, "defWeight")) fscanf(pf, "%s %lg %lg %lg %lg", tmp, &params->defWeight[0], &params->defWeight[1], &params->defWeight[2], &params->defWeight[3]);
    else if (! strcmp(field, "nThreads")) fscanf(pf, "%s %d", tmp, &params->nThreads);
    else {
       printf("Invalid input field and/ or format of input\n");
       printf("%s Doesn't match input options\n",field);
       exit(EXIT_FAILURE);
    }
    i++;
    if (i > 12) break;
  }
  fclose(pf);
  if (params->nBandStructures > 20) {printf("Too Many Band Structures: %d (up to 20 allowed)\nExiting...", params->nBandStructures); exit(EXIT_FAILURE);}
  
  if (params->nBandStructures > 1) readMultiInput(params);
  else copyBaseParamsToArrays(params);
 
  return;
 
  return;
}

/*****************************************************************************************/

void readMultiInput(param *params) {
  FILE *pf;
  int i, totKP, totAtoms, totBands;
  char field[1000], tmp[1000];

  params->nKPE = 0;
  totKP = params->nKPoints;
  totAtoms = params->nAtoms;
  totBands = params->nBands;
  if((pf = fopen("multiInput.par", "r"))==NULL){
    printf("Can't open multiInput.par\n"); exit(EXIT_FAILURE);
  }
  for (i = 0; i < params->nBandStructures; i++) {
    fscanf(pf, "%s %s %d", field, tmp, &params->na[i]);
    fscanf(pf, "%s %s %d", field, tmp, &params->ne[i]);
    fscanf(pf, "%s %s %d", field, tmp, &params->nkp[i]);
    fscanf(pf, "%s %s %lg", field, tmp, &params->mke[i]);
    params->nKPE += params->nkp[i]*params->ne[i];
    totKP -= params->nkp[i];
    totAtoms -= params->na[i];
    totBands -= params->ne[i];
  }
  fclose(pf);
  if (totBands !=0){ printf("Wrong Number of Bands! off by %d... Exiting\n", totBands); exit(EXIT_FAILURE);}
  if (totKP !=0){ printf("Wrong Number of kPoints! off by %d... Exiting\n", totKP); exit(EXIT_FAILURE);}
  if (totAtoms !=0){ printf("Wrong Number of atoms! off by %d... Exiting\n", totAtoms); exit(EXIT_FAILURE);}

  return;
}

/*****************************************************************************************/

void copyBaseParamsToArrays(param *params) {
  params->na[0] = params->nAtoms;
  params->ne[0] = params->nBands;
  params->nkp[0] = params->nKPoints;
  params->mke[0] = params->maxKE;
  params->nKPE = params->nkp[0]*params->ne[0];

  return;
}

/****************************************************************************/

void readSystem(param *params, vector *unitCellvectors, vector *defUnitCellvectors,
        atom *atoms, int fileIndex) {
  FILE *pfile;
  int i;
  char buf[1000]; 
  char tmp[1000];
  char fileName[100];
  char tmpChar[20];  
  double a,b,c;
 
  strcpy(fileName, "system_");
  sprintf(tmpChar, "%d.par", fileIndex);
  strcat(fileName, tmpChar);
  printf("opening %s to read\n",fileName);
  if ((pfile = fopen(fileName, "r") )== NULL){
    printf("Can't open %s... exiting\n", fileName);
    exit(EXIT_FAILURE);
  }

  while (fscanf(pfile, "%s", buf) == 1) {
    if (! strcmp(buf, "atoms")) {
      for (i = 0; i < params->na[fileIndex]; i++) {
        fscanf (pfile, "%s %lg %lg %lg", atoms[i].symbol, &atoms[i].pos.x, &atoms[i].pos.y, &atoms[i].pos.z);
        if (newAtomType(atoms, i)) params->nAtomTypes++; //this is done kinda weirdly... we could just imput and maybe then check
      }
    }
    else if (! strcmp(buf, "cell")) {
      for (i = 0; i < NDIM; i++) {
        fscanf(pfile, "%lg %lg %lg", &unitCellvectors[i].x, &unitCellvectors[i].y, &unitCellvectors[i].z); // unitCellvecotrs[2].z is the c parameters we need to rescale the kpoints in wz structures
      }
    }
    else if (! strcmp(buf, "scale")) fscanf(pfile, "%s %lg", tmp, &params->s[fileIndex]);
  }
  printf("Closing file\n");
  fflush(0);
  fclose(pfile);
  
  // TL: Do we need to change nElec from 2-6 to 3-5 systems? 
  params->nElec[fileIndex] = 2*params->na[fileIndex];

  // DJ: scale for deformation potential calculation
  params->defS[fileIndex] = params->s[fileIndex] + 0.001;

  /* Scale the unit cell vectors in the input file by the input parameter scale */
  for (i = 0; i < NDIM; i++) {
    unitCellvectors[i] = retScaledVector(unitCellvectors[i], params->s[fileIndex]);
    defUnitCellvectors[i] = retScaledVector(unitCellvectors[i], params->defS[fileIndex]/params->s[fileIndex]); // DJ
  }
 
  /* Scale the atom positions by the unit cell vectors */
  printf("Scaling atomic positions\n");
  for (i = 0; i < params->na[fileIndex]; i++) {
    a = atoms[i].pos.x;
    b = atoms[i].pos.y;
    c = atoms[i].pos.z;
    atoms[i].pos.x = (a * unitCellvectors[0].x + b * unitCellvectors[1].x + c * unitCellvectors[2].x);
    atoms[i].pos.y = (a * unitCellvectors[0].y + b * unitCellvectors[1].y + c * unitCellvectors[2].y);
    atoms[i].pos.z = (a * unitCellvectors[0].z + b * unitCellvectors[1].z + c * unitCellvectors[2].z);
    atoms[i].pos = retScaledVector(atoms[i].pos, 1.0);
    // DJ 
    atoms[i].defPos.x = (a * defUnitCellvectors[0].x + b * defUnitCellvectors[1].x + c * defUnitCellvectors[2].x);
    atoms[i].defPos.y = (a * defUnitCellvectors[0].y + b * defUnitCellvectors[1].y + c * defUnitCellvectors[2].y);
    atoms[i].defPos.z = (a * defUnitCellvectors[0].z + b * defUnitCellvectors[1].z + c * defUnitCellvectors[2].z);
    atoms[i].defPos = retScaledVector(atoms[i].defPos, 1.0);
  }

  params->vol[fileIndex] = getCellVolume(unitCellvectors);
  params->defVol[fileIndex] = getCellVolume(defUnitCellvectors);
  return;
}

/*****************************************************************************************/

void readExpBandStructure(double *expBandStructure, double *expDefPot, int nKPoints, int nElectrons, int fileIndex) {
  FILE *pf;
  int i, j;
  double tmp;
  char fileName[100], tmpChar[20];  
  
  strcpy(fileName, "expBandStruct_");
  sprintf(tmpChar, "%d.par", fileIndex);
  strcat(fileName, tmpChar);
  //printf("Opening %s to read\n", fileName);
  //fflush(0);
  if((pf = fopen(fileName, "r"))==NULL){ 
    printf("Cant open %s... Exiting..\n",fileName);
    exit(EXIT_FAILURE);
  }

  //printf("Reading %d k-points\n", nKPoints);
  //fflush(0);
  for (i = 0; i < nKPoints; i++) {
    //printf("k-point %d\n", i);
    //fflush(0);
    fscanf(pf, "%lg", &tmp);
    for (j = 0; j < nElectrons; j++) {
      //printf("band %d\n", i);
      //printf("%p gives %lg\n", &expBandStructure[i*nElectrons+j],expBandStructure[i*nElectrons+j] );
      //fflush(0);
      fscanf(pf, "%lg", &expBandStructure[i*nElectrons+j]); 
    }
  }
  fclose(pf);  

  // read experimental deformation potential (gap, CB, VB) in eV from file
  strcpy(fileName, "expDefPot_");
  sprintf(tmpChar, "%d.par", fileIndex);
  strcat(fileName, tmpChar);

  if ((pf = fopen(fileName, "r"))==NULL) {
    printf("Can't open %s... Exiting.\n", fileName);
    exit(EXIT_FAILURE);
  }
  fscanf(pf, "%lg %lg %lg", &expDefPot[3*fileIndex], &expDefPot[3*fileIndex+1], &expDefPot[3*fileIndex+2]);
  fclose(pf);
  
  return;
}

/****************************************************************************/

void readKPoints(vector *kPoint, vector *gVectors, int nKPoints, int fileIndex) {
  FILE *pf;
  int i;
  char fileName[100], tmpChar[20];  
  double kptx,kpty,kptz;
  kptx=0;kpty=0;kptz=0;
  strcpy(fileName, "kpoints_");
  sprintf(tmpChar, "%d.par", fileIndex);
  strcat(fileName, tmpChar);
  //printf("Opening %s to read\n", fileName);
  //fflush(0);
  if ((pf = fopen(fileName, "r"))== NULL)
  {
    printf("Can't open %s... exiting\n", fileName);
    exit(EXIT_FAILURE);

  }
  //printf("reading...\n");
  //fflush(0);
  for (i = 0; i < nKPoints; i++) {
    fscanf(pf, "%lg %lg %lg", &kptx, &kpty, &kptz);
    //kPoint[i] = retScaledVector(kPoint[i], (TWOPI / scale));
    kPoint[i].x = kptx * gVectors[0].x + kpty * gVectors[1].x + kptz * gVectors[2].x;
    kPoint[i].y = kptx * gVectors[0].y + kpty * gVectors[1].y + kptz * gVectors[2].y;
    kPoint[i].z = kptx * gVectors[0].z + kpty * gVectors[1].z + kptz * gVectors[2].z;
    kPoint[i].mag = sqrt(retDotProduct(kPoint[i], kPoint[i]));
    printf("%lg %lg %lg %lg\n",kPoint[i].x, kPoint[i].y, kPoint[i].z, kPoint[i].mag );  //for wz-structure, be careful for the z-direction!
    //fflush(0);
    }
  fclose(pf);

  return;
}

/****************************************************************************/


