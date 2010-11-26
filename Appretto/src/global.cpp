#pragma once

#include <iostream>

using namespace std;

int **glb_coord,glb_size[4],glb_vol=0;
int **loc_coord,loc_size[4],loc_vol=0;
int *glb_of_loc_ind=NULL;
int  nproc_dir[4]={0,0,0,0};
int proc_coord[4]={0,0,0,0};
int rank,rank_tot,cart_rank;

int big_endian;
const int nreals_per_spincolor=24;
const int debug=2;

MPI_Comm cart_comm;

//random generator stuff
const int ran2_ntab=32;
int *ran2_idum,*ran2_idum2,**ran2_iv,*ran2_iy;
bool random_initialized=false;

//A complex number
typedef double complex[2];

//The sum of two complex number
void complex_summ(complex a,complex b,complex c)
{
  a[0]=b[0]+c[0];
  a[1]=b[1]+c[1];
}

//Summ to the output the product of two complex number
void complex_summ_the_prod(complex a,complex b,complex c)
{
  a[0]+=b[0]*c[0]-b[1]*c[1];
  a[1]+=b[0]*c[1]+b[1]*c[0];
}

//The product of two complex number
void complex_prod(complex a,complex b,complex c)
{
  a[0]=b[0]*c[0]-b[1]*c[1];
  a[1]=b[0]*c[1]+b[1]*c[0];
}

//The product of a complex number by  the conjugate of the second
void complex_conj_prod(complex a,complex b,complex c)
{
  a[0]=b[0]*c[0]+b[1]*c[1];
  a[1]=-b[0]*c[1]+b[1]*c[0];
}

//the real amd imaginary unit
complex ONE={1,0};
complex I={0,1};

typedef complex spin[4];
typedef complex color[3];

typedef spin colorspin[3];
typedef color spincolor[4];

typedef spin spinspin[4];
typedef spinspin colorspinspin[3];

//Print a spinspin
void print_spinspin(spinspin &s)
{
  for(int id1=0;id1<4;id1++)
    {
      for(int id2=0;id2<4;id2++) cout<<s[id1][id2][0]<<","<<s[id1][id2][1]<<"\t";
      cout<<endl;
    }
}

