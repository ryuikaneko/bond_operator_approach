#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main(void)
{
  const int Lmax = 6000;
//  const int Lmax = 12000;
//  const int Lmax = 24000;
  const int Lz = 2;
  const double pi = M_PI;

  double invLmax = 1.0/Lmax;
  double invLz = 1.0/Lz;
  int i;
  int step;
//  int Nitr;
  int ix,iy;
  int minx,miny;
  double kx,ky;
  double parS,parJ,parJp,parAlpha,parJd;
  double parA,parE,parK,parLambda,parS2,parMu;
  double parR,tmp;
  double parG;
//  double mix;

  double **dispG;
  dispG = (double**)malloc(Lmax*sizeof(double*));
  dispG[0] = (double*)malloc(Lmax*Lmax*sizeof(double));
  for(i=0; i<Lmax; i++){
    dispG[i] = dispG[0] + i*Lmax;
  }

  parS = 0.5;
//  parS = 1.0;
//  parS = 1.5;
//  parS = 2.0;
  parJ = 1.0;
  parJd = 0.0;

//  Nitr = 1000;
//  mix = 0.1;

  printf("# parS J Jp Alpha Jd step  S2 Mu minGap kx/pi ky/pi  JR E(R) K(R)\n");

//for(parAlpha=0.0; parAlpha<=1.0+1.0e-10; parAlpha+=0.1){
for(parAlpha=0.0; parAlpha<=1.0+1.0e-10; parAlpha+=0.02){

//  parAlpha = 1.0;
  parJp = parAlpha*parJ;

  if(parAlpha<0.5){
    parG = 2.0 - parAlpha;
  }else{
    parG = parAlpha + 0.5/parAlpha;
  }

  for(ix=0; ix<Lmax; ix++){
    for(iy=0; iy<Lmax; iy++){
      kx = 2.0*pi*ix*invLmax;
      ky = 2.0*pi*iy*invLmax;
      dispG[ix][iy] =
        parS*(parS+1.0)*2.0/3.0*parJ*(cos(kx)+cos(ky)+parAlpha*cos(kx+ky));
    }
  }

  step = 0;

  parR = 3.0/(8.0*parS*(parS+1.0)*parJ*parG);

  parE = 0.0;
  parK = 0.0;
  for(ix=0; ix<Lmax; ix++){
    for(iy=0; iy<Lmax; iy++){
      tmp = sqrt(1.0+4.0*parR*dispG[ix][iy]);
      if(fabs(tmp)<Lmax*Lmax*Lz){
        parE += tmp;
      }
      if(fabs(1.0/tmp)<Lmax*Lmax*Lz){
        parK += 1.0/tmp;
      }
    }
  }
  parE *= invLmax*invLmax*invLz;
  parK *= invLmax*invLmax*invLz;

  parJd = parS*(parS+1.0)*parJ*parG*(20.0/3.0 - 8.0*parK);

  parS2 = 2.5 - 1.5*(parE+parK);
  parMu = - parJd*parS*(parS+1.0) + 1.5*(parE-parK)/parR;
  parA = parJd*(1-parS*(parS+1.0))-parMu;

  parLambda = sqrt(1.0+4.0*parR*dispG[Lmax/2][Lmax/2]); // gap at (pi,pi)
  minx = Lmax/2;
  miny = Lmax/2;
  for(ix=0; ix<Lmax; ix++){
    for(iy=0; iy<Lmax; iy++){
      tmp = sqrt(1.0+4.0*parR*dispG[ix][iy]);
      if(tmp<parLambda){
        parLambda = tmp;
        minx = ix;
        miny = iy;
      }
    }
  }
  parLambda *= parA/parJd;

  printf("%f %f %f %f %f %d  ",parS,parJ,parJp,parAlpha,parJd,step);
  printf("%f %f %f %f %f  ",
    parS2,parMu,parLambda,minx*2.0/Lmax,miny*2.0/Lmax);
  printf("%f %f %f\n",
    parJ*parR,parE,parK);
}

  free(dispG[0]);
  free(dispG);

  return 0;
}
