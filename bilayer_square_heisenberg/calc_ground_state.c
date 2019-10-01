#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main(void)
{
  const int Lmax = 1500;
  const int Lz = 2;
  const double pi = M_PI;

  double invLmax = 1.0/Lmax;
  double invLz = 1.0/Lz;
  int i;
  int step;
  int Nitr;
  int ix,iy;
  double kx,ky;
  double parS,parJ,parJd;
  double parA,parE,parK,parLambda,parS2,parMu;
  double parR,newR,tmpR,tmp;
  double mix;
  double conv_eps=1.0e-8;

  double **dispG;
  dispG = (double**)malloc(Lmax*sizeof(double*));
  dispG[0] = (double*)malloc(Lmax*Lmax*sizeof(double));
  for(i=0; i<Lmax; i++){
    dispG[i] = dispG[0] + i*Lmax;
  }

  parS = 0.5;
  parJ = 0.0;
  parJd = 1.0;

//  parR = 0.0;
//  newR = 0.0;

  Nitr = 1000;
  mix = 0.1;

  printf("# J Jd step  S2 Mu Gap(pi,pi)  R E(R) K(R)\n");

//for(parJd=5.0; parJd>=2.3+1.0e-10; parJd-=0.05){
for(parJ=0.0; parJ<=0.45+1.0e-10; parJ+=0.02){

  for(ix=0; ix<Lmax; ix++){
    for(iy=0; iy<Lmax; iy++){
      kx = 2.0*pi*ix*invLmax;
      ky = 2.0*pi*iy*invLmax;
      dispG[ix][iy] =
        parS*(parS+1.0)*2.0/3.0*parJ*(cos(kx)+cos(ky));
    }
  }

  parR = 0.0;
  newR = 0.0;

  for(step=0; step<Nitr; step++){
    parR = mix*newR + (1.0-mix)*parR;
    tmpR = 0.0;
    for(ix=0; ix<Lmax; ix++){
      for(iy=0; iy<Lmax; iy++){
          tmpR += 1.0/sqrt(1.0+4.0*parR*dispG[ix][iy]);
      }
    }
    newR = (2.5 - 3.0*invLmax*invLmax*invLz*tmpR)/parJd;
//    printf("%d %f %f\n",step,newR,newR-parR);
    if(fabs(parR-newR)<conv_eps) break;
    if(isnan(newR)){
      printf("# warning: newR = nan\n");
      break;
    }
  }

  parR = newR;
  parE = 0.0;
  parK = 0.0;
  for(ix=0; ix<Lmax; ix++){
    for(iy=0; iy<Lmax; iy++){
        tmp = sqrt(1.0+4.0*parR*dispG[ix][iy]);
        parE += tmp;
        parK += 1.0/tmp;
    }
  }
  parE *= invLmax*invLmax*invLz;
  parK *= invLmax*invLmax*invLz;

  parS2 = 2.5 - 1.5*(parE+parK);
  parMu = - parJd*parS*(parS+1.0) + 1.5*(parE-parK)/parR;
  parA = parJd*(1-parS*(parS+1.0))-parMu;
//  parLambda = parA*sqrt(1.0+4.0*parR*dispG[Lmax/2][Lmax/2]); // gap at (pi,pi)
  parLambda = parA/parJd*sqrt(1.0+4.0*parR*dispG[Lmax/2][Lmax/2]); // gap at (pi,pi)

  printf("%f %f %d  ",parJ,parJd,step);
  printf("%f %f %f  %f %f %f\n",
    parS2,parMu,parLambda,parR,parE,parK);
}

  free(dispG[0]);
  free(dispG);

  return 0;
}
