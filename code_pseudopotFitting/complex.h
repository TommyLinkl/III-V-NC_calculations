/****************************************************************************/
/* Include basic Libraries */
#include <stdio.h>
#include <math.h>

/****************************************************************************/
/* Define the complexnumber struct*/
typedef struct complexnumber {
  double re, im;
} complexnumber;

typedef struct cvector {
  complexnumber x, y, z;
  double mag;
} cvector;

/****************************************************************************/
/*Define some relevant functions*/
complexnumber cadd(complexnumber a, complexnumber b);
complexnumber csub(complexnumber a, complexnumber b);
complexnumber cmult(complexnumber a, complexnumber b);
complexnumber cdiv(complexnumber a, complexnumber b);
complexnumber cscale(double a, complexnumber b);
complexnumber cexpon(complexnumber a);
double cmag(complexnumber a);
complexnumber complex(double a, double b);

/*//Functions for cvectors
cvector retAddedcVectors(cvector vect1, cvector vect2);
cvector retSubtractedcVectors(cvector vect1, cvector vect2);
cvector retScaledcVector(cvector vect, complexnumber scale);
cvector retcCrossProduct(cvector vect1, cvector vect2);
cvector retZerocVector(void);
complexnumber retcDotProduct(cvector vect1, cvector vect2);
*/
