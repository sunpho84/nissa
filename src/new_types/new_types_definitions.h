#ifndef _NEW_TYPES_DEFINITIONS_H
#define _NEW_TYPES_DEFINITIONS_H

#include <mpi.h>
#include <stdint.h>

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

typedef complex color2[2];
typedef color2 su2[2];

typedef colorspinspin su3spinspin[3];

typedef complex as2t[6];
typedef su3 as2t_su3[6];

typedef MPI_Offset ILDG_Offset;
typedef MPI_File ILDG_File;

typedef momentum_t stout_pars[4];

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
//The three possibilities of quark computation
enum hmc_force_piece{GAUGE_FORCE_ONLY,QUARK_FORCE_ONLY,BOTH_FORCE_PIECES};
enum multistep_level{MACRO_STEP,MICRO_STEP};


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

struct su3_path;
struct su3_path
{
  int ivol;
  int movements;
  su3 data;
  su3_path* next;
};

//used to exponentiate for stouting
struct anti_hermitian_exp_ingredients
{
  int sign;
  double cu,su;
  double c2u,s2u;
  su3 Q,Q2;
  double u,w,theta;
  double xi0w;
  double cw;
  complex f[3];
};

//ILDG header
typedef struct
{
  uint32_t magic_no;
  uint16_t version;
  uint16_t mbme_flag;
  uint64_t data_length;
  char type[128];
} ILDG_header;

//store messages
struct ILDG_message;
struct ILDG_message
{
  bool is_last;
  char *data;
  char *name;
  ILDG_message *next;
};

//ILDG file view
typedef struct
{
  MPI_Datatype etype;
  MPI_Datatype ftype;
  MPI_Offset view_pos;
  MPI_Offset pos;
  char format[100];
} ILDG_File_view;


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
enum action_type{Wilson_action,tlSym_action};
typedef struct
{
  double beta;
  int nflavs;
  int use_bkgrd_em_field;
  quad_u1 ***backfield;
  quark_content *flav_pars;
  double E[3];
  action_type gac_type;
  double B[3];
  stout_pars stout_rho;
  int nstout_lev;
} theory_pars;

//evolution parameters for hybrid monte carlo
typedef struct
{
  double traj_length;
  double pf_action_residue;
  double md_residue;
  int nmd_steps;
  int ngauge_substeps;
  int *npseudo_fs;
} hmc_evol_pars;

typedef struct
{
  //number of hb sweeps and hits per link
  int nhb_sweeps;
  int nhb_hits;
  //the sam for overrelax
  int nov_sweeps;
  int nov_hits;
} pure_gauge_evol_pars;

typedef union
{
  hmc_evol_pars md_pars;
  pure_gauge_evol_pars pure_gauge_pars;
} evol_pars;

#endif
