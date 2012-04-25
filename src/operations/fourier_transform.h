#ifndef _FOURIER_TRANSFORM_H
#define _FOURIER_TRANSFORM_H

void Momentum(int **iP,double *bc,double *P2,double *SinP2,double **P,double **SinP,double *SinP4,int nmom);
void spincolor_FT(spincolor *S,spincolor *FT,double *theta,int **iP,int nmom);

#endif
