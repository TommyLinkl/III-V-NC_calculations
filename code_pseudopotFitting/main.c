/****************************************************************************/
/* This is the beginning of this program that fits local psudopotentials
   by comparing the resulting band structures to an expected band structure*/

/*Input files:
  input.par -- lists basic input including number of band structure to calculate
  multiInput.par -- if more than one band strcutre to be calculated list params for each
  system_X.par -- defines the crystal strcture for band structure X
  XXParams.par -- for each atom in the systems with symbol XX 
  kpoints_X.par -- gives the k-points where band structure X will be calculated (in terms of 2*pi/a)
  expBandStruct_X.par -- gives the expected band strcture for system X

Output files:
  output.dat -- gives basic summary of fitting results
  bandStruct_X.dat -- gives the energies (in eV) of the desired bands and k-points for system X
  best_XX_Y_Params.dat -- gives the best params for atom XX in system Y
  pot_XX.dat -- gives real space potential (in Hartree) corresponding to the best params for atom XX 
  basisStates.dat -- gives the plane wave basis states used for the band structure 
  iterations.dat -- gives information on each iteration of the fitting process
*/



#include "fit.h"


/****************************************************************************/
/* Initiates the main functions */

int main() {

  /****************************************************************************/
  /* Defines the variables that are used throughout this function */
  int i;
  double *energies, *expBandStructure, *expDefPot, *bandStructure, *vSpinOrbit;
  complexnumber *hamiltonian, *localPotential, *spinOrbitPot;
  vector *kPoint, *basisStates, *gVectors, *unitCellvectors, *tmpBV;
  vector *defBasisStates, *defUnitCellvectors; // DJ for deformation potential
  atom *atoms;
  param params;

  //timing
  //double tci, twi;
  
  /****************************************************************************/
  /* Allocates initial  memory */
  

  /****************************************************************************/ 
  /* Reads input parameters that determine the size of the system */

  printf("Reading Inputs\n");
  fflush(stdout);
  readInput(&params);
  int nbs = params.nBandStructures;


  printf("Allocating Memory\n");
  fflush(stdout);
  /****************************************************************************/
  /* Allocates the input dependent memory */
  if ((unitCellvectors = (vector *) calloc(nbs*NDIM, sizeof(vector)))==NULL){
    printf("error: unitCellvectors\n"); exit(EXIT_FAILURE);
  }
  if ((gVectors = (vector *) calloc(nbs*NDIM, sizeof(vector)))==NULL){
    printf("error: gVectors\n"); exit(EXIT_FAILURE);
  }
  if ((atoms = (atom *) calloc(params.nAtoms, sizeof(atom)))==NULL){
    printf("error: atoms\n"); exit(EXIT_FAILURE);
  }
  // DJ
  if ((defUnitCellvectors = (vector *) calloc(nbs*NDIM, sizeof(vector)))==NULL){
    printf("error: defUnitCellvectors\n"); exit(EXIT_FAILURE);
  }
  if ((expBandStructure = (double *) calloc(params.nKPE, sizeof(double)))==NULL){
    printf("error: expBandStructure\n"); exit(EXIT_FAILURE);
  }
  if ((expDefPot = (double *) calloc(3*params.nBandStructures, sizeof(double)))==NULL){
    printf("error: expDefPot\n"); exit(EXIT_FAILURE);
  }
  if ((bandStructure = (double *) calloc(params.nKPE, sizeof(double)))==NULL){
    printf("error: bandStructure\n"); exit(EXIT_FAILURE);
  }
  if ((kPoint = (vector *) calloc(params.nKPoints, sizeof(vector)))==NULL){
    printf("error: kPoint\n"); exit(EXIT_FAILURE);
  }
 
  /****************************************************************************/
  /* Initializes the atom positions, initial pseudopotential parameters, and calculates the basis states */
  //TODO add error handling!!! for reading inputs... I think i've mostly done this... do one last check
  printf("Reading Band Structures\n");
  fflush(stdout);

  int ina, ine, ikp, iebs;
  ina = ine = ikp = iebs = 0;

  
  for (i=0;i<nbs;i++){
    printf("Band Struct %d/%d\n",i, nbs);
    fflush(stdout);
    readExpBandStructure(&expBandStructure[iebs], expDefPot, params.nkp[i], params.ne[i], i); 
    printf("Read exp band struct\n");
    fflush(stdout);

    readSystem(&params, &unitCellvectors[i*NDIM], &defUnitCellvectors[i*NDIM], &atoms[ina], i); 
    printf("readSystem\n");
    fflush(stdout);
    //printf("ikp = %d", ikp);
    //fflush(0);
    //printf("&kPoint[ikp] = %p\n",&kPoint[ikp] );
    //printf("kPoint = %p\n",kPoint);
    getRecipSpaceBasisVectors(&unitCellvectors[i*NDIM], &gVectors[i*NDIM], &params, i);
    printf("Calculating reciprocal space basis vectors\n");
    fflush(stdout);

    readKPoints(&kPoint[ikp], &gVectors[i*NDIM], params.nkp[i], i); 
    printf("readKPoints\n");
    fflush(stdout);   


    ina += params.na[i];
    ine += params.ne[i];
    ikp += params.nkp[i];
    iebs += params.nkp[i]*params.ne[i];
    
  }

  initAtomicFormFactors(params, atoms);
  int inmbv = 0;

  printf("\nCalculating plane wave basis vectors\n");
  fflush(stdout);
  for (i=0; i<nbs; i++) {
    inmbv += cube(2*params.nmbv[i]+1);
    printf("Max number of basis vectors for system %d is %d\n", i, cube(2*params.nmbv[i]+1));
  }


  tmpBV = (vector *) calloc(inmbv, sizeof(vector));

  //TODO: the whole way we make our basis set seems a bit funky to me... 
  printf("\nFinding basis vectors under KE cuttoff\n");
  fflush(stdout);

  long inbv, inbv_squared,four_inbv_squared, ik_times_nbv_squared;
  inbv = inbv_squared = four_inbv_squared = 0;
  ik_times_nbv_squared = inmbv = 0;
  for (i = 0; i < nbs; i++) {
    calcBasisVectors(&tmpBV[inmbv], &gVectors[i*NDIM], &params, i);
    inmbv += cube(2*params.nmbv[i]+1);
    inbv += params.nbv[i];
    printf("Actual number of basis vectors for system %d is %d\n",i, params.nbv[i]);
    inbv_squared += params.nbv[i]*params.nbv[i];
    four_inbv_squared += 4*params.nbv[i]*params.nbv[i];
    ik_times_nbv_squared+=params.nbv[i]*params.nbv[i]*params.nkp[i];
  }
  printf("inbv = %ld\n", inbv);
  printf("inbv_squared = %ld\n", inbv_squared);
  printf("four_inbv_squared = %ld\n", four_inbv_squared);
  printf("ik_times_nbv_squared = %ld\n", ik_times_nbv_squared);


  printf("\nAllocating the rest of memory now we know the basis set size\n");
  fflush(stdout);
  if((vSpinOrbit = (double *) calloc(ik_times_nbv_squared,sizeof(double)))==NULL){
    printf("error: vSpinOrbit\n"); exit(EXIT_FAILURE);
  }
  if ((basisStates = (vector *) calloc(inbv, sizeof(vector)))==NULL) {
    printf("error: basisStates\n"); exit(EXIT_FAILURE);
  }
  if ((defBasisStates = (vector *) calloc(inbv, sizeof(vector)))==NULL){
    printf("error: defbasisStates\n"); exit(EXIT_FAILURE);
  }
  if ((localPotential = (complexnumber *) calloc(inbv_squared, sizeof(complexnumber)))==NULL){ //doesn't need spin
    printf("error: localPot\n"); exit(EXIT_FAILURE);
  }
  if ((hamiltonian = (complexnumber *) calloc(four_inbv_squared*params.nThreads, sizeof(complexnumber)))==NULL){
    printf("error: hamiltonian\n"); exit(EXIT_FAILURE);
  }
  if ((spinOrbitPot = (complexnumber *) calloc(four_inbv_squared*params.nThreads, sizeof(complexnumber)))==NULL){
    printf("error: spinOrbitPot\n"); exit(EXIT_FAILURE);
  }
  if ((energies = (double *) calloc(2*inbv*params.nThreads, sizeof(double)))==NULL){
    printf("error: energies\n"); exit(EXIT_FAILURE);
  }
  
  printf("Writing and sorting final basis states\n");
  fflush(stdout);  
  

  inbv = 0;
  inmbv = 0;
  for (i=0;i<nbs;i++){
    finalizeBasisStates(&basisStates[inbv], &tmpBV[inmbv], &params, i);
    inmbv += cube(2*params.nmbv[i]+1);
    inbv += params.nbv[i];
  }

  // DJ
  inbv = 0;
  for (i = 0; i < nbs; i++) {
    for (int j = 0; j < params.nbv[i]; j++) {
      defBasisStates[inbv+j] = retScaledVector(basisStates[inbv+j], params.s[i]/params.defS[i]);
    }
    inbv += params.nbv[i];
  }
  
  free(tmpBV);

  printf("Done with Basis States!\n\n");
  fflush(stdout);   

  /*****************************************************************/
  printf("Initializing spinOrbit!\n");
  fflush(stdout);
  int SOBool=0;
  for(int i =0; i<params.nAtoms;i++){
    if (atoms[i].spinOrbit != 0) {SOBool++; break;}
  }
  printf("SOBool:%i\n",SOBool );
  if (SOBool!=0){
    ina = ikp = iebs = 0;
    for (i=0;i<nbs;i++){
      calcVSpinOrbit(&vSpinOrbit[ina], &basisStates[iebs], &kPoint[ikp], params, i);
      ina+=params.nkp[i]*sqr(params.nbv[i]);
      ikp+=params.nkp[i];
      iebs+=params.nbv[i];
    }
  }
  else {
    for (int i=0;i<four_inbv_squared;i++) {spinOrbitPot[i].re=spinOrbitPot[i].im=0;}
  }
  
  /****************************************************************************/
  writeInput(atoms, params, stdout);

  printf("\nBegining Monte Carlo Run for %i iterations\n", params.nIterations);
  fflush(stdout);
  runMonteCarlo(hamiltonian, localPotential, spinOrbitPot, vSpinOrbit, basisStates, defBasisStates, expBandStructure, 
          expDefPot, kPoint, bandStructure, energies, atoms, params, SOBool, gVectors);

  /****************************************************************************/
  /* Free dynamically allocated memory */ 
  free(atoms);
  free(energies);
  free(unitCellvectors); free(gVectors);
  free(basisStates); 
  free(defUnitCellvectors); free(defBasisStates);
  free(localPotential); free(hamiltonian); 
  free(bandStructure); free(expBandStructure); free(expDefPot); free(kPoint);
  free(spinOrbitPot); free(vSpinOrbit);

  return 0;
}

/****************************************************************************/
