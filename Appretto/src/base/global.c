#pragma once

#include <stdio.h>
#include <math.h>

//nomenclature: 
//-glb is relative to the global grid
//-loc to the local one
int glb_size[4],glb_vol=0;
int loc_size[4],loc_vol=0;
//-lx is lexicografic
int **glb_coord_of_loclx;
int **loc_coord_of_loclx;
int *glblx_of_loclx=NULL;
int *glblx_of_bordlx=NULL;
int *loclx_of_bordlx=NULL;
int *dir_of_bordlx=NULL;
int *glblx_of_edgelx=NULL;
//-eor is even-odd reduced
int *loclx_parity;
int *loceo_of_loclx;
int *loclx_of_loce;
int *loclx_of_loco;
int **loce_neighup;
int **loce_neighdw;
int **loco_neighup;
int **loco_neighdw;
int loc_volr;
int appretto_eo_geom_init=0;

//neighbours of local volume + borders
int **loclx_neighdw,**loclx_neighup;

//stopping criterions for multimass inverter
int sc_standard=0,sc_unilevel=1,sc_weighted_norm2=2,sc_weighted_norm_inf=3;
int numb_known_stopping_criterion=4;
char list_known_stopping_criterion[4][1024]={"standard","unilevel","weighted_norm2","weighted_norm_inf"};

//basic mpi types
MPI_Datatype MPI_SU3;
MPI_Datatype MPI_QUAD_SU3;
MPI_Datatype MPI_SPINCOLOR;
MPI_Datatype MPI_REDSPINCOLOR;

//size of the border along the 4 dir,types for sending
int start_lx_bord_send_up[4],start_lx_bord_rece_up[4];
int start_lx_bord_send_dw[4],start_lx_bord_rece_dw[4];
int bord_dir_vol[4],bord_offset[4];
int loc_bord;
MPI_Datatype MPI_GAUGE_BORD_SEND[4],MPI_GAUGE_BORD_RECE[4];
MPI_Datatype MPI_LXSPINCOLOR_BORD_SEND[4],MPI_LXSPINCOLOR_BORD_RECE[4];
MPI_Datatype MPI_LXREDSPINCOLOR_BORD[4];

//size of the edges along the 6 directions
int edge_dir_vol[6],edge_offset[6];
int loc_edge;
int edge_numb[4][4]={{-1,0,1,2},{0,-1,3,4},{1,3,-1,5},{2,4,5,-1}};
MPI_Datatype MPI_GAUGE_EDGE_SEND[6],MPI_GAUGE_EDGE_RECE[6];

int paral_dir[4],nparal_dir;
int nproc_dir[4]={0,0,0,0};
int proc_coord[4]={0,0,0,0};
int rank_neighdw[4],rank_neighup[4];
int rank,rank_tot,cart_rank;

int big_endian;

const int nreals_per_spincolor=24;
const int nreals_per_quad_su3=72;

const int debug=1;

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
typedef color redspincolor[2];
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
void complex_subt(complex a,complex b,complex c)
{
  a[0]=b[0]-c[0];
  a[1]=b[1]-c[1];
}

//prod with real
void complex_prod_with_real(complex a,complex b,double c)
{
  a[0]=b[0]*c;
  a[1]=b[1]*c;
}

//Summ to the output the product of two complex number
//it is assumed that a!=b and a!=c

void complex_summ_the_prod(complex a,complex b,complex c)
{
  a[0]+=b[0]*c[0]-b[1]*c[1];
  a[1]+=b[0]*c[1]+b[1]*c[0];
}
void complex_subt_the_prod(complex a,complex b,complex c)
{
  a[0]-=b[0]*c[0]-b[1]*c[1];
  a[1]-=b[0]*c[1]+b[1]*c[0];
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
void complex_subt_the_conj_conj_prod(complex a,complex b,complex c)
{
  a[0]-=b[0]*c[0]-b[1]*c[1];
  a[1]-=-b[0]*c[1]-b[1]*c[0];
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
void safe_complex_conj1_prod(complex a,complex b,complex c)
{
  safe_complex_conj2_prod(a,c,b);
}

//complex prod real
void complex_prod_real(complex a,complex b,double c)
{
  a[0]=b[0]*c;
  a[1]=b[1]*c;
}

//complex prod i
void safe_complex_prod_i(complex a,complex b)
{
  double temp=b[0];
  a[0]=-b[1];
  a[1]=temp;
}
void assign_complex_prod_i(complex a){safe_complex_prod_i(a,a);}

//complex prod -i
void safe_complex_prod_minus_i(complex a,complex b)
{
  double temp=b[0];
  a[0]=b[1];
  a[1]=-temp;
}
void assign_complex_prod_minus_i(complex a){safe_complex_prod_minus_i(a,a);}

//reciprocal of a complex
void complex_reciprocal(complex rec,complex c)
{
  double module=c[0]*c[0]+c[1]*c[1];
  
  rec[0]=c[0]/module;
  rec[1]=-c[1]/module;
}

//power of a complex
void complex_pow(complex res,complex base,double exp)
{
  double module=pow(base[0]*base[0]+base[1]*base[1],exp/2);
  double anomaly=atan2(base[1],base[0])*exp;

  res[0]=module*cos(anomaly);
  res[1]=module*sin(anomaly);
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

//pi
const double PI=3.14159265358979323846;

//////////////////////////////////////////////////////////

//allocate vectors of the required length
char *allocate_vector(int length,char *tag)
{
  char *out=(char*)malloc(length);
  if(out==NULL && rank==0)
    {
      fprintf(stderr,"Error during allocation of %s\n",tag);
      MPI_Abort(MPI_COMM_WORLD,1);
    }
  return out;
}

spincolor *allocate_spincolor(int length,char *tag){return (spincolor*)allocate_vector(length*sizeof(spincolor),tag);}
redspincolor *allocate_redspincolor(int length,char *tag){return (redspincolor*)allocate_vector(length*sizeof(spincolor),tag);}
quad_su3 *allocate_quad_su3(int length,char *tag){return (quad_su3*)allocate_vector(length*sizeof(quad_su3),tag);}
su3 *allocate_su3(int length,char *tag){return (su3*)allocate_vector(length*sizeof(su3),tag);}
as2t_su3 *allocate_as2t_su3(int length,char *tag){return (as2t_su3*)allocate_vector(length*sizeof(as2t_su3),tag);}
colorspinspin *allocate_colorspinspin(int length,char *tag){return (colorspinspin*)allocate_vector(length*sizeof(colorspinspin),tag);}
su3spinspin *allocate_su3spinspin(int length,char *tag){return (su3spinspin*)allocate_vector(length*sizeof(su3spinspin),tag);}


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

//Get a spincolor from a su3spinspin
void get_spincolor_from_su3spinspin(spincolor out,su3spinspin in,int id_source,int ic_source)
{
  for(int ic_sink=0;ic_sink<3;ic_sink++)
    for(int id_sink=0;id_sink<4;id_sink++) //dirac index of sink
      {
	out[id_sink][ic_sink][0]=in[ic_sink][ic_source][id_sink][id_source][0];
	out[id_sink][ic_sink][1]=in[ic_sink][ic_source][id_sink][id_source][1];
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

//Put a spincolor into a su3spinspin
void put_spincolor_into_su3spinspin(su3spinspin out,spincolor in,int id_source,int ic_source)
{
  for(int ic_sink=0;ic_sink<3;ic_sink++)
    for(int id_sink=0;id_sink<4;id_sink++) //dirac index of sink
      {
	out[ic_sink][ic_source][id_sink][id_source][0]=in[id_sink][ic_sink][0];
	out[ic_sink][ic_source][id_sink][id_source][1]=in[id_sink][ic_sink][1];
      }
}

int min_int(int a,int b){if(a<b) return a;else return b;}
int max_int(int a,int b){if(a>b) return a;else return b;}

double min_double(double a,double b){if(a<b) return a;else return b;}
double max_double(double a,double b){if(a>b) return a;else return b;}

//swap two doubles
void swap_doubles(double *d1,double *d2)
{
  double temp=(*d1);
  (*d1)=(*d2);
  (*d2)=temp;
}

double take_time()
{
  MPI_Barrier(MPI_COMM_WORLD);
  return MPI_Wtime();
}

//Open a text file for output
FILE* open_text_file_for_output(char *outfile)
{
  FILE *fout=NULL;
  if(rank==0) fout=fopen(outfile,"w");
  if(rank==0 && fout==NULL)
    {
      fprintf(stderr,"Couldn't open the file: %s",outfile);
      MPI_Abort(MPI_COMM_WORLD,1);
    }
  return fout;
}
