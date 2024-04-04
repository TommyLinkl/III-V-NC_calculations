#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define PI     3.14159265358979323846

double aCd[4] = {-146.6990949568, 0.935852, -0.195729, 1.67983};
double aSe[4] = {22.0686801291, 4.39609, 1.20006, 0.317682};
double calcPot(double q, double *a);

int main(int argc, char const *argv[]){
	if (argc!=2) {
		printf("Invalid inputs. Usage: ./zungerPot gMag \n");
		return 1;
	}
	
	double gMag = atof(argv[1]);
	gMag*=2*PI/8.1238898;

	double vCd= calcPot(gMag, aCd);
	double vSe= calcPot(gMag, aSe);

	printf("For gMag=%.3f:\n",gMag);
	printf("vCd=%.4f\tvSe=%.4f\n",vCd,vSe);
	printf("vS=%.4f\tvA=%.4f\n",0.5*(vSe+vCd),0.5*(vCd-vSe));

	return 0;
}

double calcPot(double q, double *a) {
  double pot = (a[0]*(q*q - a[1]) / (a[2] * exp(a[3]*q*q) - 1.0));
  return pot/758.241;
}