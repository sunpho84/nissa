#pragma once

#include <stdio.h>

//nomenclature: lx is lexicografic
//glb is relative to the global grid, loc to the local one
int **glb_coord_of_loclx,glb_size[4],glb_vol=0;
int **loc_coord_of_loclx,loc_size[4],loc_vol=0;
int *glblx_of_loclx=NULL;
int *glblx_of_bordlx=NULL;
int *glblx_of_edgelx=NULL;

//neighbours of local volume + borders
int **loclx_neighdw,**loclx_neighup;

//size of the border along the 4 dir,types for sending
int bord_dir_vol[4],bord_offset[4];
int loc_bord;
MPI_Datatype MPI_SU3;
MPI_Datatype MPI_QUAD_SU3;
MPI_Datatype MPI_GAUGE_SLICE_SEND[4],MPI_GAUGE_SLICE_RECE[4];

//size of the edges along the 6 directions
int edge_dir_vol[6],edge_offset[6];
int loc_edge;
int edge_numb[4][4]={{-1,0,1,2},{0,-1,3,4},{1,3,-1,5},{2,4,5,-1}};
MPI_Datatype MPI_GAUGE_EDGE_SEND[6],MPI_GAUGE_EDGE_RECE[6];

int nproc_dir[4]={0,0,0,0};
int proc_coord[4]={0,0,0,0};
int rank_neighdw[4],rank_neighup[4];
int rank,rank_tot,cart_rank;

int big_endian;

const int nreals_per_spincolor=24;
const int nreals_per_quad_su3=72;

const int debug=3;

MPI_Comm cart_comm;

//random generator stuff
const int ran2_ntab=32;
int *ran2_idum,*ran2_idum2,**ran2_iv,*ran2_iy;
int random_initialized=0;

///////////////// New types ///////////////////

typedef double complex[2];

typedef complex spin[4];
typedef complex color[3];

typedef spin colorspin[3];
typedef color spincolor[4];

typedef spin spinspin[4];
typedef spinspin colorspinspin[3];

typedef color su3[3];
typedef su3 quad_su3[4];

typedef colorspinspin su3spinspin[3];

typedef complex as2t[6];
typedef su3 as2t_su3[6];

as2t smunu_entr[4];   //these are the sigma matrices entries
int smunu_pos[4][6];  //and positions

////////////// Operations on new types //////////////////

//The sum of two complex number
void complex_summ(complex a,complex b,complex c)
{
  a[0]=b[0]+c[0];
  a[1]=b[1]+c[1];
}

//Summ to the output the product of two complex number
//it is assumed that a!=b and a!=c
void complex_summ_the_prod(complex a,complex b,complex c)
{
  a[0]+=b[0]*c[0]-b[1]*c[1];
  a[1]+=b[0]*c[1]+b[1]*c[0];
}
void complex_summ_the_conj2_prod(complex a,complex b,complex c)
{
  a[0]+=b[0]*c[0]+b[1]*c[1];
  a[1]+=-b[0]*c[1]+b[1]*c[0];
}
void complex_summ_the_conj1_prod(complex a,complex b,complex c)
{
  complex_summ_the_conj2_prod(a,c,b);
}
void complex_summ_the_conj_conj_prod(complex a,complex b,complex c)
{
  a[0]+=b[0]*c[0]-b[1]*c[1];
  a[1]+=-b[0]*c[1]-b[1]*c[0];
}

//The product of two complex number
void unsafe_complex_prod(complex a,complex b,complex c)
{
  a[0]=b[0]*c[0]-b[1]*c[1];
  a[1]=b[0]*c[1]+b[1]*c[0];
}

//The product of a complex number by the conjugate of the second
void unsafe_complex_conj2_prod(complex a,complex b,complex c)
{
  a[0]=b[0]*c[0]+b[1]*c[1];
  a[1]=-b[0]*c[1]+b[1]*c[0];
}
void unsafe_complex_conj1_prod(complex a,complex b,complex c)
{
  unsafe_complex_conj2_prod(a,c,b);
}

//The product of the conjugate of two complex numbers
void unsafe_complex_conj_conj_prod(complex a,complex b,complex c)
{
  a[0]=b[0]*c[0]-b[1]*c[1];
  a[1]=-b[0]*c[1]-b[1]*c[0];
}

//The product of two complex number
void safe_complex_prod(complex a,complex b,complex c)
{
  double tmp=b[0]*c[0]-b[1]*c[1];
  a[1]=b[0]*c[1]+b[1]*c[0];
  a[0]=tmp;
}

//The product of a complex number by the conjugate of the second
void safe_complex_conj2_prod(complex a,complex b,complex c)
{
  double tmp=a[0]=b[0]*c[0]+b[1]*c[1];
  a[1]=-b[0]*c[1]+b[1]*c[0];
  a[0]=tmp;
}

//saturate two anti-simmetric tensors
void as2t_saturate(complex out,as2t a,as2t b)
{
  unsafe_complex_prod(out,a[0],b[0]);
  for(int munu=1;munu<6;munu++) complex_summ_the_prod(out,a[munu],b[munu]);
}

//the real amd imaginary unit
complex ONE={1,0};
complex I={0,1};

//////////////////////////////////////////////////////////

//Print a spinspin
void print_spinspin(spinspin s)
{
  for(int id1=0;id1<4;id1++)
    {
      for(int id2=0;id2<4;id2++) printf("%+016.16f,%+016.16f\t",s[id1][id2][0],s[id1][id2][1]);
      printf("\n");
    }
}

//Get a spincolor from a colorspinspin
//In a spinspin the sink index runs slower than the source                                                                                                                                                                                                                  
void get_spincolor_from_colorspinspin(spincolor out,colorspinspin in,int id_source)
{
  for(int ic_sink=0;ic_sink<3;ic_sink++)
    for(int id_sink=0;id_sink<4;id_sink++) //dirac index of sink
      {
	out[id_sink][ic_sink][0]=in[ic_sink][id_sink][id_source][0];
	out[id_sink][ic_sink][1]=in[ic_sink][id_sink][id_source][1];
      }
}

//Put a spincolor into a colorspinspin
void put_spincolor_into_colorspinspin(colorspinspin out,spincolor in,int id_source)
{
  for(int ic_sink=0;ic_sink<3;ic_sink++)
    for(int id_sink=0;id_sink<4;id_sink++) //dirac index of sink
      {
	out[ic_sink][id_sink][id_source][0]=in[id_sink][ic_sink][0];
	out[ic_sink][id_sink][id_source][1]=in[id_sink][ic_sink][1];
      }
}

int min_int(int a,int b){if(a<b) return a;else return b;}
int max_int(int a,int b){if(a>b) return a;else return b;}

double min_double(double a,double b){if(a<b) return a;else return b;}
double max_double(double a,double b){if(a>b) return a;else return b;}
