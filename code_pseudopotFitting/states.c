/****************************************************************************/
/* This files calculates the basis vectors up to a specified KE cutoff 
    given the KE cutoff and reciprocal space lattice vectors (gVectors) */

#include "fit.h"

/****************************************************************************/

void calcBasisVectors(vector *tmpBV, vector *gVectors, param *params, int fileIndex) {
  //generates nMaxBV basis vectors in the tmpBV slot and counts the ones with length under the maxKE cutoff to get nBasisVectors. 
  // &tmpBV gets filled with basis vectors
  // &params.nBasisVectors gets set to the number under the KE cutoff. 
   double MKE = params->mke[fileIndex];
   int nx, ny, nz, nbv;
   int i, j, k, index;

   nx = ny = nz = params->nmbv[fileIndex];   

   printf("GVector[0]: %f %f %f \n", gVectors[0].x, gVectors[0].y, gVectors[0].z);
   printf("GVector[1]: %f %f %f \n", gVectors[1].x, gVectors[1].y, gVectors[1].z);
   printf("GVector[2]: %f %f %f \n", gVectors[2].x, gVectors[2].y, gVectors[2].z);

   index = 0; nbv = 0;
   for (i = -nx; i < nx + 1; i++) {
      for (j = -ny; j < ny + 1; j++) {
         for (k = -nz; k < nz + 1; k++) {
            tmpBV[index] = retZeroVector();
            tmpBV[index].x += i*gVectors[0].x + j*gVectors[1].x + k*gVectors[2].x;
            tmpBV[index].y += i*gVectors[0].y + j*gVectors[1].y + k*gVectors[2].y;
            tmpBV[index].z += i*gVectors[0].z + j*gVectors[1].z + k*gVectors[2].z;
            tmpBV[index].mag = sqrt(retDotProduct(tmpBV[index], tmpBV[index]));
            if (HBAR*0.5*sqr(tmpBV[index].mag)/MASS < MKE) { 
               nbv++; //TODO: is there any way to keep hold of the indicies of the basis vectors we want ahead of time. 
            }
            index++;
         }
      }
   }
   params->nbv[fileIndex] = nbv;


   return;
}

/****************************************************************************/

void finalizeBasisStates(vector *basisStates, vector *tmpBV, param *params, int fileIndex) {
   FILE *pfile;
   int i, nbv;
   double MKE= params->mke[fileIndex];

   nbv = 0;
   for (i = 0; i < cube(2*params->nmbv[fileIndex]+1); i++) {
      if (HBAR*0.5*sqr(tmpBV[i].mag)/MASS < MKE) {
        basisStates[nbv].x = tmpBV[i].x;
        basisStates[nbv].y = tmpBV[i].y;
        basisStates[nbv].z = tmpBV[i].z;
        basisStates[nbv].mag = tmpBV[i].mag;
        nbv++;
        //TODO: should we check that we're not going over params.nBasisVectors?
      }
   }

   quickSortVectors(basisStates, 0, nbv-1);

 
   params->nqp[fileIndex] = countDistinctMagnitudes(basisStates, params->nbv[fileIndex]);

   if ((pfile = fopen("basisStates.dat", "a"))==NULL){
    printf("Error opening basisStates.dat to write... Exiting\n");
    exit(EXIT_FAILURE);
   }

   for (i = 0; i < nbv; i++) {
      fprintf(pfile, "%d % .5f % .5f % .5f % .5f\n", i, basisStates[i].x, basisStates[i].y, basisStates[i].z, basisStates[i].mag);
   }
   fclose(pfile);
   

   return;
}

/****************************************************************************/

int countDistinctMagnitudes(vector *vect, int nVectors) {
  int i, count;

  count = 1;
  for (i = 1; i < nVectors; i++) {
    if (vect[i].mag <= vect[i-1].mag) continue;
    else count++;
  }

  return count;
}

/****************************************************************************/

void quickSortVectors(vector *vectors, int l, int r) {
   int j;

   if (l < r) {
      // divide and conquer
      j = partitionVectors(vectors, l, r);
      quickSortVectors(vectors, l, j-1);
      quickSortVectors(vectors, j+1, r);
   }

   return;
}

/****************************************************************************/

int partitionVectors(vector *vectors, int l, int r) {
   int i, j;
   double pivot = vectors[l].mag;
   vector temp;

   i = l; j = r + 1;

   while(1) {
      do ++i; while (vectors[i].mag <= pivot && i <= r);
      do --j; while (vectors[j].mag > pivot);
      if (i >= j) break;
      temp = vectors[i]; vectors[i] = vectors[j]; vectors[j] = temp;
   }
   temp = vectors[l]; vectors[l] = vectors[j]; vectors[j] = temp;

   return j;
}

/****************************************************************************/
