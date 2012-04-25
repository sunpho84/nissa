#ifndef _ENDIANESS_H
#define _ENDIANESS_H

void check_endianess();
void doubles_to_doubles_changing_endianess(double *dest,double *sour,int ndoubles);
void doubles_to_floats_changing_endianess(float *dest,double *sour,int n);
void doubles_to_floats_same_endianess(float *dest,double *sour,int n);
void floats_to_doubles_changing_endianess(double *dest,float *sour,int n);
void floats_to_doubles_same_endianess(double *dest,float *sour,int n);
void floats_to_floats_changing_endianess(float *dest,float *sour,int nfloats);
void ints_to_ints_changing_endianess(int *dest,int *sour,int nints);

#endif
