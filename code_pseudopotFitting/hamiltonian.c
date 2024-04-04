/****************************************************************************/
/* this file calculates the all the matrix elements for the Hamiltonian matrix */

#include "fit.h"
//#define calcBessel(x,eps)   ((x) < (eps) ? 0.0 : (sin((x)) / ((x)*(x)) - cos((x)) / (x))) 
#define calcBessel(x,x1)   ((sin((x)) * ((x1)*(x1)) - cos((x)) * (x1))) 

/****************************************************************************/

// simply adds up the kinetic energy, local potential, and spin-orbit to get the hamiltonian matrix
void calcHamiltonianMatrix(complexnumber *hamiltonian, complexnumber *localPotential, complexnumber* spinOrbitMatrix, vector *basisStates, vector *kVector, int n, int SOBool) {
  int i, j;
  int SO_Rank = 1+SOBool;
  zeroComplexMatrix(hamiltonian, (SO_Rank)*n, (SO_Rank)*n);

  vector kPlusG;
  double preFactor;
  preFactor = ( sqr(HBAR) / (2.0 * MASS) );

  //kinetic energy
  for (i = 0; i < n; i++) {
    kPlusG = retAddedVectors(basisStates[i], *kVector);

    hamiltonian[(SO_Rank)*n*i+i].re = sqr(kPlusG.mag) * preFactor; //diagonal
    if (SOBool!=0) hamiltonian[2*n*n+n+2*n*i+i].re = sqr(kPlusG.mag) * preFactor; //diagonal
  }   
  
  // add kinetic and local potential energies to the Hamiltonian (block-diagonal)
  for (i = 0; i < n; i++) { 
    for (j = 0; j < n; j++) {
      //up up (top left) block
      hamiltonian[SO_Rank*i*n+j].re += localPotential[i*n+j].re;
      hamiltonian[SO_Rank*i*n+j].im += localPotential[i*n+j].im;
    }
  }
  if (SOBool!=0) {
    for (i = 0; i < n; i++) {
      for (j = 0; j < n; j++) {
        //dn dn (bottom right) block
        hamiltonian[2*n*n+n+2*n*i+j].re += localPotential[i*n+j].re;
        hamiltonian[2*n*n+n+2*n*i+j].im += localPotential[i*n+j].im;
      }
    }
  }
 
  if (SOBool!=0) { 
    for (i = 0; i < 4*n*n; i++) {
      hamiltonian[i].re += spinOrbitMatrix[i].re;
      hamiltonian[i].im += spinOrbitMatrix[i].im;
    }
  }

  return;
}


/****************************************************************************/

void calcPotentialEnergyMatrix(complexnumber *localPotential, vector *basisStates, atom *atoms, int nbv, int nAtoms, double cellVolume, double defCellVolume, int defBool) { 
  int i, j, k;
  double *gDiffDotTau;
  vector gDiff; 
  complexnumber *structFact; 
  //printf("Using Local PEM\n");
  int n = nbv;
  zeroComplexMatrix(localPotential, n, n);
  gDiffDotTau = (double *) calloc(nAtoms, sizeof(double));
  structFact = (complexnumber *) calloc(nAtoms, sizeof(complexnumber));
  //V_{i,j} = <G_i|V|G_j> = \sum_k [e^{-i(G_i-G_j)\cdot\tau_k} * v(|G_i-G_j|) / (V_cell)]
  printf("Cell Volume = %f\n", cellVolume);
  for (i = 0; i < n; i++) { //i is row index
    for (j = 0; j < n; j++) { //j is collumn index
      gDiff = retSubtractedVectors(basisStates[j], basisStates[i]);
      for (k = 0; k < nAtoms; k++) { 
        //TODO: Why do we divide by the number of atoms.. Its got to be wrong by just considering two atoms right on top of eachother... 
        //TODO: why do we actually store the struct factor parts? could just have temp variable, or see if needed elsewhere. not reall the structure factor either..
	      if (defBool != 0) {
          gDiffDotTau[k] = retDotProduct(gDiff, atoms[k].defPos);
          structFact[k].re = (1.0 /  ( defCellVolume )) * cos(gDiffDotTau[k]);
          structFact[k].im = (1.0 /  ( defCellVolume )) * sin(gDiffDotTau[k]); 
	      }
	      else {
          gDiffDotTau[k] = retDotProduct(gDiff, atoms[k].pos);
          structFact[k].re = (1.0 /  ( cellVolume )) * cos(gDiffDotTau[k]);
          structFact[k].im = (1.0 /  ( cellVolume )) * sin(gDiffDotTau[k]); 
	      }
        
        //structFact[k] = cscale( 1.0 / ((double)params.nAtoms * params.cellVolume), cexp(complex(0,gDiffDotTau[k])));

        // calculate pseudopotential with strain term for deformed cell
        if (defBool != 0) {
          atoms[k].fF = calcStrainPot(gDiff.mag, atoms[k].ppParams, cellVolume, defCellVolume);
        }
        else {
          atoms[k].fF = calcPot(gDiff.mag, atoms[k].ppParams);
        }
        //TODO:Does this have the wrong sign??? if inversion symmetry won't mattter?? Do we actually have inversion the way the unit cell is set up?
        //Flipped Gdiff Check this fixes
        
        // local potential energy has a delta function on spin -> block-diagonal     

        localPotential[n*i+j].re += atoms[k].fF * structFact[k].re;
        localPotential[n*i+j].im += atoms[k].fF * structFact[k].im; 
      }  
    }
  }

  free(gDiffDotTau);
  free(structFact);
  
  return;
}
/****************************************************************************/
//Calculates the reciprocial space local pseudopotential based on the given parameters at a reciprical space point with magnitude q

double calcPot(double q, double *ppParams) {
  double pot = (ppParams[0]*(q*q - ppParams[1]) / (ppParams[2] * exp(ppParams[3]*q*q) - 1.0));
  return pot;
}

/****************************************************************************/
//Calculates the reciprocial space local pseudopotential based on the given parameters at a reciprical space point with magnitude q

double calcStrainPot(double q, double *ppParams, double vol, double defVol) {
  double trEp = defVol / vol - 1.;
  // double pot = (ppParams[0]*(1. + ppParams[4]*trEp)*(q*q - ppParams[1]) / (ppParams[2] * exp(ppParams[3]*q*q) - 1.0));
  double pot = (ppParams[0]*(1. + ppParams[4]*trEp + ppParams[5]*trEp*trEp*trEp)*(q*q - ppParams[1]) / (ppParams[2] * exp(ppParams[3]*q*q) - 1.0));
  return pot;
}

/****************************************************************************/

void zeroComplexMatrix(complexnumber *mat, int n1, int n2) {
  int i, n = n1*n2;

  for (i = 0; i < n; i++) {
    mat[i].re = 0.0;
    mat[i].im = 0.0;
  }

  return;
}
/****************************************************************************/
// Calculates Vso(K,K') = 1/cellVolume * integral from 0 t0 infinity of
// dr*r^2*j1(Kr)*exp^(-(r/0.7)^2)*j1(K'r) where j1 is the 1st bessel function
// K = kpoint + basisVector and exp^(-(r/0.7)^2) is the  spin-orbit potential 
// (Hybertsen & Louie, PRB, 34, 2920 (1986))
// The nKPoint*nBasisVectors^2 of these are stored in vSpinOrbit

void calcVSpinOrbit(double *vSpinOrbit, vector *basisStates, vector *kPoint, param params, int ibs) {
  //TODO: Fix this up with OMP; check more clearly for convergence. 
  long i, j, k, gp;
  double r, eps, sum;
  long n = params.nbv[ibs];
  double dr = TWOPI/(100.0*basisStates[n-1].mag);//~=0.0089 Bohr at 25 hartree cuttoff
  //double width = 0.7; //bohr
  double  width = params.SOWidth;
  double rCut = sqrt(width*width*16.0*log(10.0));//~=4.2488 bohr; V(rCut)= 1e-16
  int nCut = (int)(rCut/dr); //around 400 maybe?
  vector gikp, gjkp;


  eps = dr / 10.0;
  //printf("Can use %d Threads\n", omp_get_max_threads());
  omp_set_num_threads(params.nThreads);
  #pragma omp parallel for private(k,i,j,gikp,gjkp,sum,gp,r)
  for (k = 0; k < params.nkp[ibs]; k++) {
    //if(omp_get_thread_num()==0){printf("Using %d Threads\n", omp_get_num_threads()); fflush(0);}
    //printf("Using %d Threads (thread %d)\n", omp_get_num_threads(), omp_get_thread_num()); fflush(0);
    //printf("Starting k=%ld (thread %d)\n", k, omp_get_thread_num()); fflush(0);
    for (i = 0; i < n; i++) {
      for (j = 0; j <= i; j++) {
        gikp = retAddedVectors(basisStates[i], kPoint[k]);
        gjkp = retAddedVectors(basisStates[j], kPoint[k]);
        sum = 0.0;
        // perform numerical integration 
        for (gp = 1; gp < nCut; gp++) {
          r = ((double)(gp) * dr);
          sum += (sqr(r) * calcBessel(gikp.mag*r, (1.00/(gikp.mag*r + EPS)) ) 
            * exp(-(sqr(r/width))) * calcBessel(gjkp.mag*r, (1.00/(gjkp.mag*r + EPS))) * dr);                       
        }
        vSpinOrbit[k*n*n+i*n+j] = sum / (params.vol[ibs]);
        vSpinOrbit[k*n*n+j*n+i] = sum / (params.vol[ibs]);
      }
    }
    //printf("done with k=%ld (thread %d)\n", k, omp_get_thread_num()); fflush(0);
  }

  return;
}
/****************************************************************************/
void calcSpinOrbitPotentialMatrix(complexnumber *spinOrbitMatrix, double *vSpinOrbit, vector *basisStates, vector kVector, int kindex, atom *atoms, param params,int ibs) {
  int i, j, k;
  double *gDiffDotTau, preFactor, *lambda;
  vector ibsPlusK, jbsPlusK, gCrossProduct, gDiff;
  vector upUpRe, upUpIm, upDownRe, upDownIm, downUpRe, downUpIm, downDownRe, downDownIm;
  complexnumber gcpDotUpUp, gcpDotDownDown, gcpDotDownUp, gcpDotUpDown;
  complexnumber *structFact;

  gDiffDotTau = (double *) calloc(params.na[ibs], sizeof(double));
  structFact = (complexnumber *) calloc(params.na[ibs], sizeof(complexnumber));
  lambda = (double *) calloc(params.na[ibs], sizeof(double));

  // -i*<s|SpinMatrices|s'> vectors: //Should I represent these as 8 3-vectors? four for real and for for imag
  /*
  sMat= (vector *) calloc(8,sizeof(vector));
  for(i=0; i<8;i++){
    
  }
  */
  upUpRe = retZeroVector();
  upUpIm.x = 0.0, upUpIm.y = 0.0, upUpIm.z = -0.5;
  downDownRe = retZeroVector();
  downDownIm.x = 0.0, downDownIm.y = 0.0, downDownIm.z = 0.5;
  upDownRe.x = 0.0, upDownRe.y = -0.5, upDownRe.z = 0.0;
  upDownIm.x = -0.5, upDownIm.y = 0.0, upDownIm.z = 0.0;
  downUpRe.x = 0.0, downUpRe.y = 0.5, downUpRe.z = 0.0;
  downUpIm.x = -0.5, downUpIm.y = 0.0, downUpIm.z = 0.0;
  

  int n = params.nbv[ibs];
  zeroComplexMatrix(spinOrbitMatrix, 2*n, 2*n); 
  
  for (i = 0; i < n; i++) { // i = k
    for (j = 0; j < n; j++) { // j = k'
      gDiff = retSubtractedVectors(basisStates[j], basisStates[i]);
      // TODO (Mar 22 2019): check if gDiff should be j-i instead of i-j
      ibsPlusK = retAddedVectors(basisStates[i], kVector);
      jbsPlusK = retAddedVectors(basisStates[j], kVector);    
      // avoid dividing by 0 in spinOrbitPrefactor ... should be fine b/c vSpinOrbit will be zero
      if (ibsPlusK.mag < EPS || jbsPlusK.mag < EPS) {
        continue; //This should be OK because the V_SO=0 for gDiff=0
      }
      else {
        preFactor = 12.0 * PIE / (ibsPlusK.mag * jbsPlusK.mag);
        gCrossProduct = retCrossProduct(ibsPlusK, jbsPlusK);
        gcpDotUpUp.re = retDotProduct(upUpRe, gCrossProduct); // is always 0
        gcpDotUpUp.im = retDotProduct(upUpIm, gCrossProduct);
        gcpDotDownDown.re = retDotProduct(downDownRe, gCrossProduct); // is always 0
        gcpDotDownDown.im = retDotProduct(downDownIm, gCrossProduct);
        gcpDotUpDown.re = retDotProduct(upDownRe, gCrossProduct);
        gcpDotUpDown.im = retDotProduct(upDownIm, gCrossProduct);
        gcpDotDownUp.re = retDotProduct(downUpRe, gCrossProduct);
        gcpDotDownUp.im = retDotProduct(downUpIm, gCrossProduct); 
        for (k = 0; k < params.na[ibs]; k++) { 
          gDiffDotTau[k] = retDotProduct(gDiff, atoms[k].pos);
          structFact[k].re =  cos(gDiffDotTau[k]);
          structFact[k].im =  sin(gDiffDotTau[k]);
          lambda[k] = atoms[k].spinOrbit * vSpinOrbit[kindex*n*n+i*n+j];
         
          //TODO: tidy up, check if we need to do this for each atom or if we can do the structure factor as a single number...
          //Add the structure factor for each atom to the matrix
          spinOrbitMatrix[2*n*i+j] = cadd(cscale(preFactor * lambda[k], cmult(gcpDotUpUp, structFact[k])), spinOrbitMatrix[2*n*i+j]);
          spinOrbitMatrix[2*n*n+n+2*n*i+j] = cadd(cscale(preFactor * lambda[k], cmult(gcpDotDownDown, structFact[k])), spinOrbitMatrix[2*n*n+n+2*n*i+j]);
          spinOrbitMatrix[n+2*n*i+j] = cadd(cscale(preFactor * lambda[k], cmult(gcpDotUpDown, structFact[k])), spinOrbitMatrix[n+2*n*i+j]);
          spinOrbitMatrix[2*n*n+2*n*i+j] = cadd(cscale(preFactor * lambda[k], cmult(gcpDotDownUp, structFact[k])), spinOrbitMatrix[2*n*n+2*n*i+j]);
        }
        // TODO: check if dot product in outside or inside of loop over the atoms matters
      }
    }
  }

  free(gDiffDotTau);
  free(structFact);
  free(lambda);

  return;
}
