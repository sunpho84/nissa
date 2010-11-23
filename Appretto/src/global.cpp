#pragma once

int **glb_coord,glb_size[4],glb_vol=0;
int **loc_coord,loc_size[4],loc_vol=0;
int *glb_of_loc_ind=NULL;
int  nproc_dir[4]={0,0,0,0};
int proc_coord[4]={0,0,0,0};
int rank,rank_tot,cart_rank;

int big_endian;
const int debug=1;

MPI_Comm cart_comm;

//random generator stuff
const int ran2_ntab=32;
int *ran2_idum,*ran2_idum2,**ran2_iv,*ran2_iy;
bool random_initialized=false;

typedef double complex[2];

typedef complex spin[4];
typedef complex color[3];

typedef spin colorspin[3];
typedef color spincolor[4];

typedef spin spinspin[4];
typedef spinspin colorspinspin[3];
