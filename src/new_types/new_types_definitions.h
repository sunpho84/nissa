#ifndef _NEW_TYPES_DEFINITIONS_H
#define _NEW_TYPES_DEFINITIONS_H

#include <mpi.h>
#include <stdint.h>

#include <lemon.h>

#include "../base/macros.h"

///////////////// New types ///////////////////

typedef int coords[4];
typedef double momentum_t[4];

typedef int intpair[2];
typedef uint32_t checksum[2];

typedef double complex[2];
typedef complex quad_u1[4];

typedef complex spin[4];
typedef complex color[3];

typedef spin colorspin[3];
typedef color spincolor[4];

typedef colorspin spincolorspin[4];
typedef spincolorspin colorspincolorspin[3];

typedef spin spinspin[4];
typedef spinspin colorspinspin[3];

typedef color su3[3];
typedef su3 quad_su3[4];

typedef colorspinspin su3spinspin[3];

typedef complex as2t[6];
typedef su3 as2t_su3[6];

//this is just for avoid misleading, but is nothing more that a spinspin
typedef complex spin1field[4];
typedef spin1field spin1prop[4];

//quadruple precision float
typedef double float_128[2];
typedef float_128 complex_128[2];
typedef complex_128 color_128[3];
typedef color_128 spincolor_128[4];

//Random types
enum rnd_type{RND_ALL_PLUS_ONE,RND_ALL_MINUS_ONE,RND_UNIF,RND_Z2,RND_Z4,RND_GAUSS};
//Source type
enum source_type{POINT_SOURCE,UNDILUTED_SOURCE,COLOR_DILUTED_SOURCE,SPIN_DILUTED_SOURCE,SPINCOLOR_DILUTED_SOURCE};

///////////////////// New structures //////////////

//The structure for the random generator
typedef struct
{
  int idum;
  int idum2;
  int iv[ran2_ntab];
  int iy;
} rnd_gen;

//The structure for gamma matrix
typedef struct
{
  int pos[4];
  complex entr[4];
} dirac_matr;

//nissa vector
struct nissa_vect;

struct nissa_vect
{
  int nel;
  int size_per_el;
  
  char tag[nissa_vect_string_length];
  char type[nissa_vect_string_length];
  
  char file[nissa_vect_string_length];
  int line;
  
  nissa_vect *prev;
  nissa_vect *next;
  
  uint32_t flag;
  
  //padding to keep memory alignment
  char pad[nissa_vect_alignment-(3*sizeof(int)+2*sizeof(void*)+3*nissa_vect_string_length+sizeof(uint32_t))%nissa_vect_alignment];
};

//nissa file reader
typedef struct
{
  int open;
  int reading;
  int nbytes_per_site;
  char *buf;
  MPI_File *reader_file;
  LemonReader *lemon_reader;
} nissa_reader;

//rational approximation
typedef struct
{
  char name[20];
  double minimum;
  double maximum;
  double cons;
  double exp_power;
  int nterms;
  double *poles;
  double *weights;
} rat_approx;

//quark content
typedef struct
{
  int deg;
  double mass;
  double re_pot;
  double im_pot;
  double charge;
} quark_content;

//theory content
typedef struct
{
  double beta;
  int nflavs;
  quad_u1 ***backfield;
  quark_content *flav_pars;
} theory_pars;

//evolution parameters
typedef struct
{
  double traj_length;
  double pf_action_residue;
  double md_residue;
  int nmd_steps;
  int ngauge_substeps;
  int *npseudo_fs;
} evol_pars;

#endif
