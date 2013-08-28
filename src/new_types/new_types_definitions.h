#ifndef _NEW_TYPES_DEFINITIONS_H
#define _NEW_TYPES_DEFINITIONS_H

#ifdef HAVE_CONFIG_H
 #include "config.h"
#endif

#include <mpi.h>
#include <math.h>
#include <stdint.h>
#include <map>
#include <vector>
#include <string>

#include "../base/macros.h"

#ifdef SPI
 #include <spi/include/kernel/MU.h>
#endif

///////////////// New types ///////////////////

typedef int coords[4];
typedef double momentum_t[4];

typedef int intpair[2];
typedef uint32_t checksum[2];

typedef double complex[2];
typedef complex quad_u1[4];

typedef complex spin[4];
typedef complex color[3];
typedef complex halfspin[2];

typedef spin colorspin[3];
typedef color spincolor[4];
typedef color halfspincolor[2];

typedef colorspin spincolorspin[4];
typedef spincolorspin colorspincolorspin[3];

typedef spin spinspin[4];
typedef spinspin colorspinspin[3];

typedef color su3[3];
typedef su3 quad_su3[4];
typedef su3 oct_su3[8];

typedef complex color2[2];
typedef color2 su2[2];

typedef colorspinspin su3spinspin[3];

typedef complex as2t[6];
typedef su3 as2t_su3[6];

typedef su3 squared_staples_t[4][6];

typedef squared_staples_t rectangular_staples_t;

#ifdef BGQ
 typedef complex bi_complex[2];
 typedef bi_complex bi_color[3];
 typedef bi_color bi_su3[3];
 typedef bi_su3 bi_oct_su3[8];
 typedef bi_color bi_spincolor[4];
 typedef bi_color bi_halfspincolor[2];
 typedef bi_complex bi_halfspin[2];

 #ifdef BGQ_EMU
  typedef bi_complex reg_bi_complex;
 #else
  typedef vector4double reg_bi_complex;
 #endif
#endif

#ifdef USE_MPI_IO
typedef MPI_Offset ILDG_Offset;
typedef MPI_File ILDG_File;
#else
typedef off_t ILDG_Offset;
typedef FILE* ILDG_File;
#endif

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

typedef uint8_t coords_5D[5];

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
#include "../base/vectors.h"

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

////////////////////// optimized 2 pts //////////////////////

//match of forward and backward components id
struct forw_back_comp_id_t: std::pair<uint8_t,uint8_t>
{
  uint8_t &forw() {return this->first;}
  uint8_t &back() {return this->second;}
  forw_back_comp_id_t(uint8_t forw_comp_id,uint8_t back_comp_id)
    {
      forw()=forw_comp_id;
      back()=back_comp_id;
    }
};

//list of correlation function to which the combo contributes, and weight
struct corr_id_weight_t: std::map<uint16_t,double> {};

//Structure useful to compute a list of correlation function of the type
// [Re/Im]Tr[G_sink*S_forw*G_sour*G5*S_back^+*G5]
//we order S as sink-source-re/im.
//To compute each correlation function we need to summ 32 contribution,
//relative to the product of the real and imaginary components.
//Two main patterns are possible, according to the the fact that we must 
//compute the real or imaginary part of the correlation function,
//and if one of the two dirac matrices are imaginary.
//We compact the product of the various correlation functions, so that
//we can compute only once the repeated products.
struct two_pts_comp_t: std::map<forw_back_comp_id_t,corr_id_weight_t>
{
  void add_sink_source_corr(uint16_t corr_id,double weight,int re_im,uint8_t sink_igamma,uint8_t sour_igamma);
  void print(FILE *fout);
  void summ_the_loc_forw_back_contractions(double *out,double *S_forw,double *S_back,int nel,int twall);
  void scream();
  std::vector<int> pattern_list;
};

//////////////////////////////////////////////////////////////////////

// Type to hold the position of linear algebra output data (see "two_stage_computations" doc for explenations)
struct two_stage_computation_pos_t
{
  int *inter_fr_in_pos; //offset for intermediate result
  int *final_fr_inter_pos; //offset for final result from intermediate
  int *inter_fr_recv_pos; //offset for intermediate from nissa_recv_buf
  two_stage_computation_pos_t()
  {
    inter_fr_in_pos=final_fr_inter_pos=inter_fr_recv_pos=NULL;
  }
  void free()
  {
    if(inter_fr_in_pos!=NULL) nissa_free(inter_fr_in_pos);
    if(final_fr_inter_pos!=NULL) nissa_free(final_fr_inter_pos);
    if(inter_fr_recv_pos!=NULL) nissa_free(inter_fr_recv_pos);
  }
  ~two_stage_computation_pos_t() {free();}
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

//parameters to measure all rectangles path
struct all_rect_meas_pars_t
{
  int flag;
  char path[1024];
  
  //hyp pars
  double hyp_alpha0,hyp_alpha1,hyp_alpha2;
  
  //ape pars
  double ape_alpha;
  int nape_levls,*nape_iters;
  
  //intervals for rectangles
  int Tmin,Tmax,Dmin,Dmax;
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

//parameters to ape smear
struct ape_pars_t
{
  int nlev;
  double alpha;
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

//The structure to hold trajectory statistics
struct hmc_traj_stat_el_t
{
  int n;
  double sx;
  double s2x;
hmc_traj_stat_el_t(int n,double sx,double s2x) : n(n),sx(sx),s2x(s2x) {}
hmc_traj_stat_el_t() : n(0),sx(0),s2x(0) {}
  hmc_traj_stat_el_t operator+=(double x){n++;sx+=x;s2x+=x*x;return (*this);}
  void ave_err(double &ave,double &err){ave=sx;err=0;if(n>=2){ave=sx/n;err=s2x/n;err=sqrt((err-ave*ave)/(n-1));}}
};
typedef std::pair<int,int> hmc_traj_stat_pars_t;
typedef std::map<hmc_traj_stat_pars_t,hmc_traj_stat_el_t> hmc_traj_stat_t;

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

//results of a unitarity check
struct unitarity_check_result_t
{
  int nbroken_links;
  double average_diff;
  double max_diff;
};

//parameters to adapt the algorithm
struct adaptative_algorithm_pars_t
{
  int use_for;
  hmc_traj_stat_t stat;
  hmc_traj_stat_pars_t current;
  void init_from_text(char *text);
  std::string save_to_text();
  void add(int a,int b,double x) {stat[std::make_pair(a,b)]+=x;}
  void set(int a,int b,int n,double sx,double s2x) {stat[std::make_pair(a,b)]=hmc_traj_stat_el_t(n,sx,s2x);}
  void set_current(hmc_evol_pars_t &a) {current=std::make_pair(a.nmd_steps,a.ngauge_substeps);}
};

//parameters for pure gauge theory
struct pure_gauge_evol_pars_t
{
  //number of hb sweeps and hits per link
  int nhb_sweeps;
  int nhb_hits;
  //the same for overrelax
  int nov_sweeps;
  int nov_hits;
};

union evol_pars_t
{
  hmc_evol_pars_t hmc_evol_pars;
  pure_gauge_evol_pars_t pure_gauge_evol_pars;
};

//out and in buffer
struct comm_t
{
  //bgq specific structures, in alternative to ordinary MPI
#ifdef SPI
  //descriptors
  MUHWI_Descriptor_t *descriptors;
  MUHWI_Destination spi_dest[8];
#else
  //destinations and source ranks
  int send_rank[8],recv_rank[8];
  //requests and message
  MPI_Request requests[16];
  int nrequest,imessage;
#endif
  
  //communication in progress
  int comm_in_prog;
  //local size
  size_t nbytes_per_site;
  //size of the message
  uint64_t tot_mess_size;
  //offsets
  int send_offset[8],message_length[8],recv_offset[8];
};

#endif
