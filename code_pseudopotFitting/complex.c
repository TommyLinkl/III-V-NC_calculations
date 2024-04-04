#include "complex.h"
#include <math.h>

#define complexI (complexnumber) {0.00,1.00};

complexnumber cadd(complexnumber a, complexnumber b){
	complexnumber sum;
	sum.re=a.re+b.re;
	sum.im=a.im+b.im;
	return sum;
}
complexnumber csub(complexnumber a, complexnumber b){
	complexnumber dif;
	dif.re=a.re-b.re;
	dif.im=a.im-b.im;
	return dif;
}
complexnumber cmult(complexnumber a, complexnumber b){
	complexnumber prod;
	prod.re=(a.re * b.re) - (a.im * b.im);
	prod.im=(a.re * b.im) + (a.im * b.re);
	return prod;
}
complexnumber cdiv(complexnumber a, complexnumber b){
	complexnumber quot;
	quot.re=(a.re*b.re + a.im*b.im)/ (b.re*b.re + b.im*b.im);
	quot.im=(a.im*b.re - a.re*b.im)/ (b.re*b.re + b.im*b.im);
	return quot;
}
complexnumber cscale(double a, complexnumber b){
	complexnumber out;
	out.re = a * b.re;
	out.im = a * b.im;
	return out;
}

complexnumber cexpon(complexnumber a){
	complexnumber out;
	out.re = exp(a.re)*cos(a.im);
	out.im = exp(a.re)*sin(a.im);
	return out;
}

double cmag(complexnumber a){
	double mag = sqrt(a.re*a.re+a.im*a.im);
	return mag;
}

complexnumber complex(double a, double b){
	complexnumber out;
	out.re = a;
	out.im = b;
	return out;
}

//cVector Functions
/*
cvector retAddedcVectors(cvector vect1, cvector vect2){
	cvector vect;
	vect.x = cadd(vect1.x,vect2.x); vect.y = cadd(vect1.y,vect2.y); vect.z = cadd(vect1.z,vect2.z);
  	vect.mag = sqrt(cmag(retcDotProduct(vect, vect)));
  	return vect;
}

cvector retSubtractedcVectors(cvector vect1, cvector vect2){
	cvector vect;
	vect.x = csub(vect1.x,vect2.x); vect.y = csub(vect1.y,vect2.y); vect.z = csub(vect1.z,vect2.z);
  	vect.mag = sqrt(cmag(retcDotProduct(vect, vect)));
  	return vect;
}

cvector retScaledcVector(cvector vectin, complexnumber scale){
	cvector vect;
	vect.x = cscale(scale,vectin.x); vect.y = cscale(scale,vectin.y); vect.z = cscale(scale,vectin.z);
  	vect.mag = cscale(vectin.mag);
  	return vect;
}
cvector retcCrossProduct(cvector vect1, cvector vect2){
	cvector crossProduct;

  	crossProduct.x = csub(cmult(vect1.y , vect2.z) , cmult(vect1.z , vect2.y));
  	crossProduct.y = csub(cmult(vect1.z , vect2.x) , cmult(vect1.x , vect2.z));
  	crossProduct.z = csub(cmult(vect1.x , vect2.y) , cmult(vect1.y , vect2.x));
  	crossProduct.mag = sqrt(cmag(retcDotProduct(crossProduct, crossProduct)));

  return crossProduct;
}
cvector retZerocVector(void){
	cvector zerovec;
	zerovec.x.re = zerovec.x.im = 0;
	zerovec.y.re = zerovec.y.im = 0;
	zerovec.y.re = zerovec.y.im = 0;
	zerovec.mag = 0;
}
complexnumber retcDotProduct(cvector vect1, cvector vect2){
	return csum( csum( cmult(vect1.x,vect2,x) , cmult(vect1.y,vect2.y) ) cmult(vect1.z,vect2.z));
}
*/