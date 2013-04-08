#ifndef _NEW_TYPES_DEFINITIONS_H
#define _NEW_TYPES_DEFINITIONS_H

#ifdef HAVE_CONFIG_H
 #include "config.h"
#endif

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

typedef su3 squared_staples_t[4][6];

typedef squared_staples_t rectangular_staples_t;

typedef MPI_Offset ILDG_Offset;
typedef MPI_File ILDG_File;

//this is just for avoid misleading, but is nothing more that a spinspin
typedef complex spin1field[4];
typedef spin1field spin1prop[4];

//quadruple precision float
typedef double float_128[2];
typedef float_128 complex_128[2];
typedef complex_128 color_128[3];
typedef color_128 spincolor_128[4];

//octpuple
typedef double float_256[4];
typedef double float_256_unr[5];

//Random types
enum rnd_t{RND_ALL_PLUS_ONE,RND_ALL_MINUS_ONE,RND_UNIF,RND_Z2,RND_Z4,RND_GAUSS};
//Source type
enum source_t{POINT_SOURCE,UNDILUTED_SOURCE,COLOR_DILUTED_SOURCE,SPIN_DILUTED_SOURCE,SPINCOLOR_DILUTED_SOURCE};
//The three possibilities of quark computation
enum hmc_force_piece{GAUGE_FORCE_ONLY,QUARK_FORCE_ONLY,BOTH_FORCE_PIECES};
enum multistep_level{MACRO_STEP,MICRO_STEP};

///////////////////// New structures ////////////////////

//The structure for the random generator
struct rnd_gen
{
  int idum;
  int idum2;
  int iv[ran2_ntab];
  int iy;
};

//The structure for gamma matrix
struct dirac_matr
{
  int pos[4];
  complex entr[4];
};

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
  double c0,c1;
  double cu,su;
  double c2u,s2u;
  su3 Q,Q2;
  double u,w,theta;
  double xi0w;
  double cw;
  complex f[3];
};

//ILDG header
struct ILDG_header
{
  uint32_t magic_no;
  uint16_t version;
  uint16_t mbme_flag;
  uint64_t data_length;
  char type[128];
};

//store messages
struct ILDG_message
{
  bool is_last;
  char *data;
  char *name;
  ILDG_message *next;
};

//ILDG file view
struct ILDG_File_view
{
  MPI_Datatype etype;
  MPI_Datatype ftype;
  MPI_Offset view_pos;
  MPI_Offset pos;
  char format[100];
};

//rational approximation
struct rat_approx_t
{
  char name[20];
  double minimum;
  double maximum;
  double cons;
  int num,den;
  double exp_power;
  int degree;
  double *poles;
  double *weights;
};

//quark content
struct quark_content_t
{
  int deg;
  double mass;
  double re_pot;
  double im_pot;
  double charge;
};

//parameters to compute the chiral condensate
struct chiral_cond_pars_t
{
  int flag;
  char path[1024];
  double residue;
  int nhits;
};

//parameters to compute the magnetization
struct magnetization_pars_t
{
  int flag;
  char path[1024];
  double residue;
  int nhits;
};

//parameters to compute time pseduoscalar correlator
struct pseudo_corr_pars_t
{
  int flag;
  char path[1024];
  double residue;
  int nhits;
};

//parameters to measure topology properties
struct top_meas_pars_t
{
  int flag;
  char path[1024];
  int cool_nsteps;
  int cool_overrelax_flag;
  double cool_overrelax_exp;
  int meas_each;
};

typedef momentum_t stout_coeff_t[4];

//structure to store data for stouting
struct stout_link_staples
{
  su3 C;
  su3 Omega;
  su3 Q;
};

//parameters to stout
struct stout_pars_t
{
  int nlev;
  stout_coeff_t rho;
};

//background em field parameters
struct em_field_pars_t
{
  int flag;
  double E[3];
  double B[3];
};

//theory content
enum gauge_action_name_t{Wilson_action,tlSym_action};
struct theory_pars_t
{
  double beta;
  int nflavs;
  quad_u1 ***backfield;
  quark_content_t *quark_content;
  gauge_action_name_t gauge_action_name;
  stout_pars_t stout_pars;
  em_field_pars_t em_field_pars;
  chiral_cond_pars_t chiral_cond_pars;
  magnetization_pars_t magnetization_pars;
  pseudo_corr_pars_t pseudo_corr_pars;
};

//evolution parameters for hybrid monte carlo
struct hmc_evol_pars_t
{
  int skip_mtest_ntraj;
  double traj_length;
  double pf_action_residue;
  double md_residue;
  int nmd_steps;
  int ngauge_substeps;
  int *npseudo_fs;
};

//parameters for pure gauge theory
struct pure_gauge_evol_pars_t
{
  //number of hb sweeps and hits per link
  int nhb_sweeps;
  int nhb_hits;
  //the sam for overrelax
  int nov_sweeps;
  int nov_hits;
};

union evol_pars_t
{
  hmc_evol_pars_t hmc_evol_pars;
  pure_gauge_evol_pars_t pure_gauge_evol_pars;
};

//////////////////////////////////////// BGQ specifics ///////////////////////////////////

#ifdef BGQ

#include <spi/include/kernel/MU.h>

//////////////// new types /////////////////

//type to hold the 5D coordinates
typedef uint8_t coords_5D[5];

//structure used to hold spi buffers 
struct spi_comm_t
{
  //communication in progress
  int comm_in_prog;
  //size of the buffers, buffers
  uint64_t buf_size;
  char *send_buf,*recv_buf;
  //counter for received bytes
  volatile uint64_t recv_counter;
  //descriptors
  MUHWI_Descriptor_t *descriptors;
  //bat
  MUSPI_BaseAddressTableSubGroup_t spi_bat_gr;
  uint32_t bat_id[2];
};

#endif

#endif
