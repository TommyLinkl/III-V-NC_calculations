/****************************************************************************/
/* This files calculates the band structure along specified directions in k-space */
//TODO: Add acceptance ratio, 

#include "fit.h"

/****************************************************************************/

void runMonteCarlo(complexnumber *hamiltonian, complexnumber *potential, complexnumber *spinOrbitPot, double *vSpinOrbit, vector *basisStates, 
        vector *defBasisStates, double *expBandStructure, double *expDefPot, vector *kPoint, double *bandStructure, double *energies, atom *atoms,
        param params,int SOBool, vector *gVectors) {

  FILE *pf;
  int i;
  double bestMSE, currentMSE, newMSE;
  double *tmpPseudos; /* Used to hold current pseudos while testing new ones */
  double *bestPseudos; 
  double defPot[3*params.nBandStructures]; // gap, CB, VB
  //timing
  double tci, twi;

  tmpPseudos = (double *) calloc(NPParams*params.nAtoms, sizeof(double));
  bestPseudos = (double *) calloc(NPParams*params.nAtoms, sizeof(double));
 
  pf = fopen("iterations.dat", "w");
  //writeCurrentTime(pf);
  //writeSeparation(pf); 

  tci = (double)clock(); 
  twi = (double)time(NULL); 

  calcBandStructure(bandStructure, hamiltonian, potential, spinOrbitPot, vSpinOrbit, basisStates, kPoint, energies, atoms, params, SOBool);  
  calcDeformationPotentials(defPot, hamiltonian, potential, spinOrbitPot, basisStates, defBasisStates, energies, atoms, params);

  bestMSE = currentMSE = calcAllIndividualMSE(bandStructure, expBandStructure, kPoint, defPot, expDefPot, atoms, params, gVectors);
  savePseudoParams(tmpPseudos, atoms, params); savePseudoParams(bestPseudos, atoms, params);
  writeIteration(currentMSE, 0, atoms, params, defPot, pf);
  fflush(0);

  for (i = 0; i < params.nIterations; i++) {
    if ((i+1) % 50 == 0) printf("Iteration %d \n", i);
    //TODO: I've hard coded some approximate scaling of step to param scale. Could be optimized better
    updatePseudoParams(atoms, params); 
   
    calcBandStructure(bandStructure, hamiltonian, potential, spinOrbitPot, vSpinOrbit, basisStates, kPoint, energies, atoms, params, SOBool);
    calcDeformationPotentials(defPot, hamiltonian, potential, spinOrbitPot, basisStates, defBasisStates, energies, atoms, params);

    //TODO: I've added a weighting of bands near the middle of desired range. Not very cleanly done and could break if params.ne[i] are not all 8 
    //remover after updated
    // static int first  = 0;
    // if (first ==0) printf("You still havent updated/fixed wieghting function in fit.c:48");
    // first++;
    newMSE = calcAllIndividualMSE(bandStructure, expBandStructure, kPoint, defPot, expDefPot, atoms, params, gVectors); 
    
    if ((i+1) % 1 == 0) writeIteration(newMSE, i+1, atoms, params, defPot, pf);
    if (newMSE < bestMSE) { 
      bestMSE = currentMSE = newMSE;
      savePseudoParams(bestPseudos, atoms, params); savePseudoParams(tmpPseudos, atoms, params);
    }
    else if (exp(-params.beta * (sqrt(newMSE) - sqrt(currentMSE))) > (0.5*(randNumber()+1.0))) {
       currentMSE = newMSE; 
       savePseudoParams(tmpPseudos, atoms, params);
    } 
    else { revertPseudoParams(atoms, tmpPseudos, params); } 
  }
  if (params.nIterations!=0){
    revertPseudoParams(atoms, bestPseudos, params);
    calcBandStructure(bandStructure, hamiltonian, potential, spinOrbitPot, vSpinOrbit, basisStates, kPoint, energies, atoms, params, SOBool);    
    bestMSE = calcAllIndividualMSE(bandStructure, expBandStructure, kPoint, defPot, expDefPot, atoms, params, gVectors);
  }
  printf("Done fitting, writing output\n");
  fflush(0);
  writeResults(bandStructure, kPoint, bestMSE, atoms, params);
  writeCurrentTime(pf);
  fclose(pf);

  // Calculate and print the final r-space local pseudopotentials
  calcRealSpacePotential(atoms, params);

  // Free dynamically allocated memory
  free(tmpPseudos);
  free(bestPseudos);

  return;
}

/****************************************************************************/

void calcBandStructure(double *bandStructure, complexnumber *hamiltonian, complexnumber *potential, complexnumber *spinOrbitPot, double *vSpinOrbit,
                  vector *basisStates, vector *kPoint, double *energies, atom *atoms, param params, int SOBool) {
  
  int i, j, k;
  int inbv = 0;
  int inbv_squared = 0;
  int four_inbv_squared = 0;
  int ikp = 0;
  int ina = 0;
  int ikpe = 0;
  long inkp_sqr_BV = 0;
  double tci_comp[4] = {0,0,0,0};
  double twi_comp[4] = {0,0,0,0};
  double tci,twi;

  for (i = 0; i < params.nBandStructures; i++) {

    calcPotentialEnergyMatrix(&potential[inbv_squared], &basisStates[inbv], &atoms[ina], params.nbv[i], params.na[i], params.vol[i], params.defVol[i], 0);

    //parallelize over k-points
    omp_set_num_threads(params.nThreads); 
    omp_set_max_active_levels(1);

    #pragma omp parallel for private(j,k, tci,twi,tci_comp,twi_comp)

    for (j = 0; j < params.nkp[i]; j++) {
      long tid = omp_get_thread_num();
      zeroComplexMatrix(&hamiltonian[four_inbv_squared*params.nThreads+4*sqr(params.nbv[i])*tid], 2*params.nbv[i], 2*params.nbv[i]);

      if (SOBool!=0) calcSpinOrbitPotentialMatrix(&hamiltonian[four_inbv_squared*params.nThreads+4*sqr(params.nbv[i])*tid], &vSpinOrbit[inkp_sqr_BV], &basisStates[inbv], kPoint[j+ikp], j, &atoms[ina], params, i);
      
      calcHamiltonianMatrix(&hamiltonian[four_inbv_squared*params.nThreads+4*sqr(params.nbv[i])*tid], &potential[inbv_squared], &spinOrbitPot[four_inbv_squared], &basisStates[inbv], &kPoint[j+ikp], params.nbv[i], SOBool);

      diagonalizeHamiltonian(&hamiltonian[four_inbv_squared*params.nThreads+4*sqr(params.nbv[i])*tid], &energies[2*inbv*params.nThreads+2*params.nbv[i]*tid], (1+SOBool) * params.nbv[i]);
    
      if (SOBool==0) for (k = params.nbv[i]; k>0; k--) energies[2*inbv*params.nThreads+2*params.nbv[i]*tid+2*k-1]=energies[2*inbv*params.nThreads+2*params.nbv[i]*tid+2*k-2] = energies[2*inbv*params.nThreads+2*params.nbv[i]*tid+k-1];
          
      for (k = 0; k < params.ne[i]; k++) bandStructure[ikpe + j*params.ne[i] + k] = energies[2*inbv*params.nThreads+2*params.nbv[i]*tid+k]*AUTOEV;
    }

    inbv += params.nbv[i];
    inbv_squared += params.nbv[i]*params.nbv[i];
    four_inbv_squared += 4*params.nbv[i]*params.nbv[i];
    ikp += params.nkp[i];
    ina += params.na[i];
    ikpe += params.nkp[i]*params.ne[i];
    inkp_sqr_BV+=params.nkp[i]*sqr(params.nbv[i]);
  }

  return;
}

/****************************************************************************/

void calcBandGaps(double *cBands, double *vBands, complexnumber *hamiltonian, 
        complexnumber *potential, complexnumber *spinOrbitPot, vector *basisStates, 
        double *energies, atom *atoms, param params, double *cellVolumes, double *defCellVolumes, int defBool) {
  // Note that here cBands[] and vBands[] are in units of Hartree, same as energies[]. 
  int inbv = 0;
  int inbv_squared = 0;
  int four_inbv_squared = 0;
  int ina = 0;
  vector gammaPoint[1];

  gammaPoint[0] = retZeroVector();

  for (int i = 0; i < params.nBandStructures; i++) {

    calcPotentialEnergyMatrix(&potential[inbv_squared], &basisStates[inbv], &atoms[ina], params.nbv[i], params.na[i], cellVolumes[i], defCellVolumes[i], defBool);
    calcHamiltonianMatrix(&hamiltonian[four_inbv_squared], &potential[inbv_squared], &spinOrbitPot[four_inbv_squared],
            &basisStates[inbv], &gammaPoint[0], params.nbv[i], 0);
    diagonalizeHamiltonian(&hamiltonian[four_inbv_squared], energies, 1*params.nbv[i]);
    cBands[i] = energies[params.nElec[i]];
    vBands[i] = energies[params.nElec[i]-1];
    
    inbv += params.nbv[i];
    inbv_squared += params.nbv[i]*params.nbv[i];
    four_inbv_squared += 4*params.nbv[i]*params.nbv[i];
    ina += params.na[i];

  }

  return;

}

/****************************************************************************/

void calcDeformationPotentials(double *defPot, 
        complexnumber *hamiltonian, complexnumber *potential, complexnumber *spinOrbitPot,
        vector *basisStates, vector *defBasisStates, double *energies, atom *atoms, param params) {

  double cBands[params.nBandStructures], defCBands[params.nBandStructures];
  double vBands[params.nBandStructures], defVBands[params.nBandStructures];

  calcBandGaps(cBands, vBands, hamiltonian, potential, spinOrbitPot, basisStates, energies, atoms, params, params.vol, params.defVol, 0);
  calcBandGaps(defCBands, defVBands, hamiltonian, potential, spinOrbitPot, defBasisStates, energies, atoms, params, params.vol, params.defVol, 1);

  for (int i = 0; i < params.nBandStructures; i++) {
    defPot[3*i+1] = (cBands[i]-defCBands[i])/(params.vol[i]-params.defVol[i])*0.5*(params.vol[i]+params.defVol[i]) * AUTOEV; // CB
    defPot[3*i+2] = (vBands[i]-defVBands[i])/(params.vol[i]-params.defVol[i])*0.5*(params.vol[i]+params.defVol[i]) * AUTOEV; // VB
    defPot[3*i] = defPot[3*i+1]-defPot[3*i+2]; // gap
    // printf("%i %.8f %.8f\n", i, params.vol[i], params.defVol[i]);
    // printf("%i %.8f %.8f\n", i, bandGaps[i], defBandGaps[i]);
    // printf("%i gap, CB, VB: %.8f %.8f %.8f\n", i, (cBands[i]-vBands[i])*AUTOEV, cBands[i]*AUTOEV, vBands[i]*AUTOEV);
    // printf("%i defGap, defCB, defVB: %.8f %.8f %.8f\n", i, (defCBands[i]-defVBands[i])*AUTOEV, defCBands[i]*AUTOEV, defVBands[i]*AUTOEV);
    // printf("%i defPot, defPotCB, defPotVB: %.8f %.8f %.8f\n", i, defPot[3*i]*AUTOEV, defPot[3*i+1]*AUTOEV, defPot[3*i+2]*AUTOEV); 
  }

  return;

}
    
/****************************************************************************/

void updatePseudoParams(atom *atoms, param params) {
  int i, j;
  //Scale the step size to be more in line with the sensitivity of each param
  //double paramStepSizeScale[NPParams];
  // paramStepSizeScale[0] = 0;// 30.00; //a0
  // paramStepSizeScale[1] = 0;// 1;  //a1
  // paramStepSizeScale[2] = 0;// 1;  //a2
  // paramStepSizeScale[3] = 0;// 0.1;  //a3
  // paramStepSizeScale[4] = 1; // a4

  for (i = 0; i < params.nAtoms; i++) {
    if (newAtomType(atoms, i)) { //TODO why is this inside the for loop??
      if (! strcmp(atoms[i].symbol, "Ga")) { // if Ga, MC step on all 6 params 
	      for (j = 0; j < NPParams; j++) {atoms[i].ppParams[j] += params.stepSizeScale[j]*params.stepSize*randNumber();} /* randNumber returns a number between -1.0 and 1.0 */
      }
      else { // if not Ga, MC step only on a0-a4, leaving a5 intact as 0
	      for (j = 0; j < NPParams-1; j++) {atoms[i].ppParams[j] += params.stepSizeScale[j]*params.stepSize*randNumber();} /* randNumber returns a number between -1.0 and 1.0 */
      }
    } 
    else { 
      copyPreviousAtomParams(atoms, i); 
    }
  }

  return;
}

/****************************************************************************/

void copyPreviousAtomParams(atom *atoms, int i) {
  int j, k;

  for (j = 0; j < i; j++) {
    if (! strcmp(atoms[j].symbol, atoms[i].symbol)) {
      for (k = 0; k < NPParams; k++) {
        atoms[i].ppParams[k] = atoms[j].ppParams[k];
      }
      return;
    }
  }
  printf("Error copying atoms params for atom %s: %d\n Exiting...\n", atoms[i].symbol, i );
  exit(EXIT_FAILURE);
  return; // should never reach here
}

/****************************************************************************/

void revertPseudoParams(atom *atoms, double *tmpPseudos, param params) {
  int i, j;

  for (i = 0; i < params.nAtoms; i++) {
    for (j = 0; j < NPParams; j++) atoms[i].ppParams[j] = tmpPseudos[i*NPParams+j];
  }

  return;
}

/****************************************************************************/

void savePseudoParams(double *tmpPseudos, atom *atoms, param params) {
  int i, j;

  for (i = 0; i < params.nAtoms; i++) {
    for (j = 0; j < NPParams; j++)  tmpPseudos[i*NPParams+j] = atoms[i].ppParams[j];
  }

  return;
}

/****************************************************************************/

double calcAllIndividualWeightedMSE(double *bandStructure, double *expBandStructure, param params) {
  int i, j; 
  double mse[params.nBandStructures], totalMSE = 0.0;
  double bandWeight[params.nBands];
  printf("MSE\n");
  fflush(0);
  int ikpe = 0;
  for (i = 0; i < params.nBandStructures; i++) {
    // Weight bands closest to band-edge more than those far apart
    for (j = 0; j < params.ne[i]/2; j++) {
      if (j < params.ne[i]/8) bandWeight[j] = bandWeight[params.ne[i]-j-1] = 0.06125; 
      else if (j < params.ne[i]/4) bandWeight[j] = bandWeight[params.ne[i]-j-1] = 2.0;
      else bandWeight[j] = bandWeight[params.ne[i]-j-1] = 4.0;
    }
    // Calculated weighted mean squared error between bandStructure and expBandStructure
    mse[i] = calcWeightedMSE(&bandStructure[ikpe], &expBandStructure[ikpe], &bandWeight[0], params.ne[i], params.nkp[i]);
    totalMSE += mse[i];
    ikpe += params.ne[i]*params.nkp[i];
  }
  
  return (totalMSE / (double)(params.nBandStructures));
}

/*****************************************************************************************/

double calcWeightedMSE(double *array1, double *array2, double *weights, int nElectrons, int nKPoints) {
  int i, j, index, n = nElectrons*nKPoints; 
  double mse = 0.0;
  
  for (i = 0; i < nKPoints; i++) {
    for (j = 0; j < nElectrons; j++) {
      index = i*nElectrons+j; 
      if ((array2[index] > EPS) || (array2[index] < -EPS)) {
        mse += weights[j] * sqr((array1[index] - array2[index]));
      }
      else {
        n--;
      }
    }
  }
  
  return (mse / (double)(n));
}

/*****************************************************************************************/

double calcAllIndividualMSE(double *bandStructure, double *expBandStructure, vector *kPoint, double *defPot, double *expDefPot, atom *atoms, param params, vector *gVectors) {
  int i,j; 
  double mse[params.nBandStructures], totalMSE = 0.0, defMSE = 0.0;
  int ikp = 0;
  int ikpe = 0;

  for (i = 0; i < params.nBandStructures; i++) {
    for (j = 0; j < params.nkp[i]; j++){
      // mse[i] = calcMSE(&bandStructure[ikpe], &expBandStructure[ikpe], &kPoint[ikp+j], params.ne[i]);
      mse[i] = calcMSE(&bandStructure[ikpe], &expBandStructure[ikpe], &kPoint[ikp+j], params.ne[i], &gVectors[i*NDIM]);
      totalMSE += mse[i];
      ikpe += params.ne[i];
    }
    ikp += params.nkp[i];
    /***
    // Add weighted contribution from deformation potential
    defMSE += params.defWeight[0]*sqr(defPot[3*i] - expDefPot[3*i]);
    defMSE += params.defWeight[1]*sqr(defPot[3*i+1] - expDefPot[3*i+1]);
    defMSE += params.defWeight[2]*sqr(defPot[3*i+2] - expDefPot[3*i+2]);
    ***/
    defMSE += params.defWeight[0]*sqr(defPot[3*i] - expDefPot[3*i]) / sqr(expDefPot[3*i]);
    defMSE += params.defWeight[1]*sqr(defPot[3*i+1] - expDefPot[3*i+1]) / sqr(expDefPot[3*i+1]);
    defMSE += params.defWeight[2]*sqr(defPot[3*i+2] - expDefPot[3*i+2]) / sqr(expDefPot[3*i+2]);
    
  }
  // Add penalty for large a4 values
  for (i = 0; i < params.nAtoms; i++) {
      if (newAtomType(atoms, i)) {
          defMSE += params.defWeight[3]*sqr(atoms[i].ppParams[4]);
      }
  }

  return (defMSE + totalMSE / (double)(params.nKPoints));
}

/*****************************************************************************************/

double calcMSE(double *array1, double *array2, vector *kVector, int arrayLength, vector *gVectors) {
  int i, count; 
  count = 0;
  double mse = 0.0;
  double kweight, bandWeight;
  vector L_point, G_point, X_point; 
  double Lptx,Lpty,Lptz,Gptx,Gpty,Gptz,Xptx,Xpty,Xptz;
  //kweight = exp(-2.*sqr((*kVector).mag));
  
  Lptx = 0.5; Lpty = 0.5; Lptz = 0.5;
  Gptx = 0.0; Gpty = 0.0; Gptz = 0.0; 
  Xptx = 0.5; Xpty = 0.0; Xptz = 0.5; 

  L_point.x = Lptx * gVectors[0].x + Lpty * gVectors[1].x + Lptz * gVectors[2].x;
  L_point.y = Lptx * gVectors[0].y + Lpty * gVectors[1].y + Lptz * gVectors[2].y;
  L_point.z = Lptx * gVectors[0].z + Lpty * gVectors[1].z + Lptz * gVectors[2].z;
  L_point.mag = sqrt(retDotProduct(L_point, L_point));
  G_point.x = Gptx * gVectors[0].x + Gpty * gVectors[1].x + Gptz * gVectors[2].x;
  G_point.y = Gptx * gVectors[0].y + Gpty * gVectors[1].y + Gptz * gVectors[2].y;
  G_point.z = Gptx * gVectors[0].z + Gpty * gVectors[1].z + Gptz * gVectors[2].z;
  G_point.mag = sqrt(retDotProduct(G_point, G_point));
  X_point.x = Xptx * gVectors[0].x + Xpty * gVectors[1].x + Xptz * gVectors[2].x;
  X_point.y = Xptx * gVectors[0].y + Xpty * gVectors[1].y + Xptz * gVectors[2].y;
  X_point.z = Xptx * gVectors[0].z + Xpty * gVectors[1].z + Xptz * gVectors[2].z;
  X_point.mag = sqrt(retDotProduct(X_point, X_point));
  
  //printf("Ax %.10f, Ay %.10f, Az %.10f, Mx %.10f, My %.10f, Mz %.10f, ", (A_point.x),(A_point.y),(A_point.z), (M_point.x),(M_point.y),(M_point.z));
  kweight = exp(-10.*sqr((*kVector).mag)) + exp(-10.*sqr((retSubtractedVectors(*kVector, L_point)).mag)) + exp(-10.*sqr((retSubtractedVectors(*kVector, G_point)).mag))  + exp(-10.*sqr((retSubtractedVectors(*kVector, X_point)).mag));
  if ((*kVector).mag<=EPS || retSubtractedVectors(*kVector, L_point).mag<=EPS || retSubtractedVectors(*kVector, G_point).mag<=EPS || retSubtractedVectors(*kVector, X_point).mag<=EPS) {kweight=500.0;} 
  
  //printf("kx %.6f, ky %.6f, kz %.6f, k-A %.6f, k-M %.6f, kweight %.4f\n", (*kVector).x, (*kVector).y, (*kVector).z , retSubtractedVectors(*kVector, A_point).mag, retSubtractedVectors(*kVector, M_point).mag ,kweight);
  
  //kweight = 1.0;
  //printf("kmag %.8f, kweight %.8f\n", (*kVector).mag, kweight);
  for (i = 0; i < arrayLength; i++) {
    if (i<2) bandWeight = 0.0; 
    else if (i==2 || i==3) bandWeight = 40.0; 
    else if (i>3 & i<8) bandWeight = 100.0; 
    else if (i==8) bandWeight = 350.0; 
    else if (i==9) bandWeight = 180.0; 
    else if (i==10) bandWeight = 40.0; 
    else if (i==11) bandWeight = 30.0; 
    else if (i==12) bandWeight = 10.0; 
    else if (i>12) bandWeight = 0.0; 
    
    if ((array2[i] > EPS) || (array2[i] < -EPS)) {
      //mse += (i) * (arrayLength-1-i) * sqr((array1[i] - array2[i]));
      mse += bandWeight * sqr((array1[i] - array2[i]));
      count++;
    }
  }
  
  if (count == 0) return 0;

  return (kweight * mse / (double)(count));
}
/****************************************************************************/

void storeEnergies(double *bandStructure, double *energies, int index, int length) {
  int j;

  for (j = 0; j < length; j++) {
    bandStructure[index*length+j] = energies[j]*AUTOEV;
  }

  return;
}

/****************************************************************************/
