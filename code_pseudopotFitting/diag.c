/****************************************************************************/
/* This file diagnolizes the hamiltonian matrix and stores the eigenvalues in energies */

#include "fit.h"


/****************************************************************************/

void diagonalizeHamiltonian(complexnumber *hamiltonian, double *energies, int nBV) {
  long n = nBV;
  //int lwork = 3*n;
  int info, layout;
  char  jobz, uplo;

  MKL_Complex16* a = (MKL_Complex16*) &hamiltonian[0];

  /* Compute only eigenvalues */
  jobz = 'N';
  uplo = 'U';
  info = 0;
  layout  = LAPACK_ROW_MAJOR;

  omp_set_num_threads(1);
  info = LAPACKE_zheev(layout, jobz, uplo, n, a,n,&energies[0]);
  
  if (info) { // 0 gets returned on successful exit
    printf("Error in zheev_ in diag.c. Paramater %d is incorrect!\n", -info);
    printf("Program is exiting!!!\n");
    exit(EXIT_FAILURE);
  }

  return;
}

/****************************************************************************/

void rapidDiagonalizeHamiltonian(complexnumber *hamiltonian, double *energies){
  
  int layout  = LAPACK_ROW_MAJOR;
  char jobz = 'N';
  char range = 'I';
  char uplo = 'U';
  int n = 1;
  MKL_Complex16* a = (MKL_Complex16*) &hamiltonian[0];
  int lda = n;
  double vl = 0;
  double vu = 0; //unused for range = 'I'
  int il = 1;
  int iu = 1;
  double abstol = 1e-10;
  lapack_int* m = 0;
  double* w = &energies[0];
  MKL_Complex16* z = (MKL_Complex16*) &hamiltonian[0];
  int ldz = 1;
  lapack_int* isuppz = (lapack_int*) calloc(n, sizeof(lapack_int));
  
  printf("begining diag...\n");
  fflush(0);
  int info = LAPACKE_zheevr(layout, jobz, range, uplo, n, a, lda, vl, vu, il, iu, abstol, m, w, z, ldz, isuppz);

  if (info) { // 0 gets returned on successful exit
    printf("Error in zheev_ in diag.c. Paramater %d is incorrect!\n", -info);
    printf("Program is exiting!!!\n");
    exit(EXIT_FAILURE);
  }
  printf("Found %lld eigs\n", *m);
  fflush(0);

}