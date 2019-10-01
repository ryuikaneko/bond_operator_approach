#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main(void)
{
//  const int Lmax = 2500;
  const int Lmax = 1500;
//  const int Lmax = 100;
  const int Lz = 2;
  const double pi = M_PI;

  double invLmax = 1.0/Lmax;
  double invLz = 1.0/Lz;
  int i;
  int step;
  int Nitr;
  int ix,iy;
  double kx,ky;
  double parS,parJ,parJp,parJd;
  double parS2,parMu,parA;
  double newS2,newMu;
  double tmpS2,tmpMu;
  double mix;

  double **dispE;
  dispE = (double**)malloc(Lmax*sizeof(double*));
  dispE[0] = (double*)malloc(Lmax*Lmax*sizeof(double));
  for(i=0; i<Lmax; i++){
    dispE[i] = dispE[0] + i*Lmax;
  }
  double **dispLambda;
  dispLambda = (double**)malloc(Lmax*sizeof(double*));
  dispLambda[0] = (double*)malloc(Lmax*Lmax*sizeof(double));
  for(i=0; i<Lmax; i++){
    dispLambda[i] = dispLambda[0] + i*Lmax;
  }

/*
  parS = 0.5;
  parJ = 1.0;
  parJp = 0.0;
  parJd = 10.0;
//  parJd = 3.0;
//  parJd = 2.287;
//  parJd = 1.0;
*/

  parS = 0.5;
//  parJ = 0.0;
//  parJ = 0.2;
  parJ = 0.4;
//  parJ = 0.4375;
  parJp = 0.0;
  parJd = 1.0;

  parS2 = -0.5;
  parMu = -0.75;
  newS2 = -0.5;
  newMu = -0.75;

  Nitr = 1000;
//  mix = 0.05;
  mix = 0.1;
//  mix = 0.2;
//  mix = 0.5;
//  mix = 0.8;

//for(parJd=2.0; parJd<=2.5+1.0e-10; parJd+=0.1){
  for(step=0; step<Nitr; step++){

    parS2 = mix*newS2 + (1.0-mix)*parS2;
    parMu = mix*newMu + (1.0-mix)*parMu;

    parA = parJd*(1-parS*(parS+1.0))-parMu;
    for(ix=0; ix<Lmax; ix++){
      for(iy=0; iy<Lmax; iy++){
        kx = 2.0*pi*ix*invLmax;
        ky = 2.0*pi*iy*invLmax;
        dispE[ix][iy] =
          parS*(parS+1.0)*parS2*2.0/3.0*(
          + parJ*( cos(kx)+cos(ky) )
          + parJp*( cos(kx+ky) )
          );
        dispLambda[ix][iy] =
          sqrt(fabs(parA*(parA+4.0*dispE[ix][iy])));
      }
    }

    tmpS2 = 0.0;
    tmpMu = 0.0;
    for(ix=0; ix<Lmax; ix++){
      for(iy=0; iy<Lmax; iy++){
        if(dispLambda[ix][iy] > 1.0e-10){
          tmpS2 += (parA+2.0*dispE[ix][iy])/dispLambda[ix][iy];
          tmpMu += dispE[ix][iy]/dispLambda[ix][iy];
        }
      }
    }
//    newS2 = 2.5 - 3.0*invLmax*invLmax*tmpS2;
//    newMu = - parJd*parS*(parS+1.0) + 6.0*parA*invLmax*invLmax/parS2*tmpMu;
    newS2 = 2.5 - 3.0*invLmax*invLmax*invLz*tmpS2;
    newMu = - parJd*parS*(parS+1.0) + 6.0*parA*invLmax*invLmax*invLz/parS2*tmpMu;

    printf("%d %f %f %f: %f %f\n",
      step,newS2,newMu,dispLambda[Lmax/2][Lmax/2],
      parS2-newS2,parMu-newMu);

    if(fabs(parS2-newS2)<1.0e-12 && fabs(parMu-newMu)<1.0e-12) break;

  }
    printf("%f ",parJd);
    printf("%d %f %f %f\n",
      step,newS2,newMu,dispLambda[Lmax/2][Lmax/2]);
//}

  free(dispE[0]);
  free(dispE);
  free(dispLambda[0]);
  free(dispLambda);

  return 0;
}
