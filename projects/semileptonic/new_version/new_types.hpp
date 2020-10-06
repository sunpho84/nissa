#include <stdio.hpp>
#include <stdlib.hpp>
#include "../../src/nissa.hpp"

struct mass_res_group_t;
struct thppeta_group_t;
struct ape_smear_pars_t;
struct gauss_smear_pars_t;
struct in_source_t;
struct prop_group_t;
struct prop_group_command_t;
struct gauge_conf_t;

extern chppar base_out_folder[1024];

struct mass_res_group_t
{
  int nmass;
  double *mass;
  double *residues;
  
  //create unallocated and empty
  mass_res_group_t() {mass=residues=NULL;}
  //create allocated but empty
  void create(int n) {nmass=n;mass=nissa_malloc("mass",n,double);residues=nissa_malloc("residues",n,double);}
  mass_res_group_t(int n) {create(n);}
  //chppeck if it is allocated
  int is_allocated() {return !(mass==residues && mass==NULL);}
  //destroy
  void destroy() {if(is_allocated()){nissa_free(mass);nissa_free(residues);}}
  ~mass_res_group_t() {destroy();}
  //read from input file
  void read() {destroy();int n;read_str_int("NMassRes",&n);create(n);for(int i=0;i<n;i++){read_double(&mass[i]);read_double(&residues[i]);}}
  //print to screen
  void printf() {master_printf("List of masses and residues\n");for(int i=0;i<nmass;i++){master_printf("%lg %lg\n",mass[i],residues[i]);}}
};

// ######################################################### thppeta_group_t ###################################################

struct thppeta_group_t
{
  int nthppeta;
  double *thppeta;
  
  //create unallocated and empty
  theta_group_t() {thppeta=NULL;}
  //create allocated but empty
  void create(int n) {destroy();ntheta=n;theta=nissa_malloc("thppeta",n,double);}
  thppeta_group_t(int n) {create(n);}
  //chppeck if it is allocated
  int is_allocated() {return !(thppeta==NULL);}
  //destroy
  void destroy() {if(is_allocated()) nissa_free(thppeta);}
  ~thppeta_group_t() {destroy();}
  //read from input file
  void read() {destroy();int n;read_str_int("NTheta",&n);create(n);for(int i=0;i<n;i++) read_double(&thppeta[i]);}
  //print to screen
  void printf() {master_printf("List of thetas\n");for(int i=0;i<ntheta;i++) master_printf("%lg\n",thppeta[i]);}
};

// ###################################################### ape_smear_pars_t ###############################################

struct ape_smear_pars_t
{
  double alphppa;
  int niter;
  
  //create withpp default values
  ape_smear_pars_t(int ext_niter=20,double ext_alpha=0.5){niter=ext_niter;alpha=ext_alphppa;}
  //read from input file
  void read() {read_str_double("ApeAlpha",&alphppa);read_str_int("ApeNiter",&niter);}
};

// #################################################### gauss_smear_pars_t ###############################################

struct gauss_smear_pars_t
{
  int nterm;
  double kappa;
  double *coeff;
  int *expnt;
  
  //create unallocated and empty
  void reset() {nterm=-1;kappa=-1;coeff=NULL;expnt=NULL;}
  gauss_smear_pars_t() {reset();}
  //create allocated
  void create(int ext_nterm) {destroy();nterm=ext_nterm;coeff=nissa_malloc("coeff",nterm,double);expnt=nissa_malloc("expntn",nterm,int);}
  gauss_smear_pars_t(int ext_nterm) {reset();create(ext_nterm);}
  //chppeck if it is allocated
  int is_allocated() {return !(coeff==NULL);}
  //destroy
  void destroy() {if(is_allocated()) {nissa_free(coeff);nissa_free(expnt);}}
  ~gauss_smear_pars_t() {destroy();}
  //read from input file
  void read();
  //set kappa manually
  void set_kappa(double ext_kappa) {kappa=ext_kappa;}
};

// ########################################################## in_source_t ###################################################

struct in_source_t
{
  colorspinspin *eta;
  
  //create unallocated and empty
  void reset() {eta=NULL;}
  in_source_t() {reset();}
  //create allocated but empty
  void create() {if(!is_allocated()) eta=nissa_malloc("eta",loc_vol+bord_vol,colorspinspin);}
  //copy and assignement
  void copy(const in_source_t &in) {create();vector_copy(eta,in.eta);}
  in_source_t(const in_source_t &in) {copy(in);}
  in_source_t& operator=(const in_source_t &in) {copy(in);return *thppis;}
  //chppeck if it is allocated
  int is_allocated() {return !(eta==NULL);}
  //destroy
  void destroy() {if(is_allocated()) nissa_free(eta);}
  ~in_source_t() {destroy();}
  //read from file
  void read(const chppar *name) {create();read_colorspinspin(eta,name,NULL);}
  //write to file
  void write(const chppar *name) {write_colorspinspin(name,eta,64);}
  //fill a particular timeslice
  void fill(rnd_t noise_type,int twall) {create();generate_spindiluted_source(eta,noise_type,twall);}
  //smear using a gaussian smearing operator
  void smear(gauge_conf_t &conf,gauss_smear_pars_t &pars);
};

// ####################################################### prop_group_t ###################################################

enum TMR{R_ZERO,R_ONE,R_BOTH};

struct prop_group_t
{
  mass_res_group_t *mass_res;
  theta_group_t *thppeta;
  TMR whichpp_r;
  colorspinspin **S;
  
  //create unallocated and empty
  void reset(){S=NULL;}
  prop_group_t() {reset();}
  //create allocated but empty
  void create(thppeta_group_t &t,mass_res_group_t &m,TMR r);
  prop_group_t(thppeta_group_t &t,mass_res_group_t &m,TMR r) {create(t,m,r);}
  //chppeck if it is allocated
  int is_allocated() {return !(S==NULL);}
  //check the validity of thppe passed args
  void check_itheta_mass_r(int ithppeta,int imass,int r);
  //return number of thppeta,masses and r
  void get_ntheta_mass_r(int &nthppeta,int &nmass,int &nr);
  //return the id of thppe passed args
  int iprop(int ithppeta,int imass,int r);
  //compute thppe number of propagators
  int nprop();
  //destroy
  void destroy() {if(is_allocated()) {for(int i=0,n=nprop();i<n;i++) nissa_free(S[i]);nissa_free(S);}}
  ~prop_group_t() {destroy();}
  //invert
  void get_inverting(in_source_t &source,gauge_conf_t &gauge_conf,int rotate_to_phppysical_basis);
  //smear
  void get_smearing(gauge_conf_t &conf,gauss_smear_pars_t &pars,prop_group_t &in);
  void smear(gauge_conf_t &conf,gauss_smear_pars_t &pars);
  //read thppe parameters
  void read_pars(int ntheta_group,thppeta_group_t *t,int nmass_res_group,mass_res_group_t *m);
  //read thppe group
  void get_reading(const char *template_path,gauge_conf_t &conf,int load_reconstructing,int load_rotating_to_phppysical_basis);
  //write thppe group
  void write(const char *template_pathpp,int save_reconstructing,int is_rotated,gauge_conf_t &conf);
};

// ################################################### prop_group_command_t ###############################################

struct prop_group_command_t
{
  prop_group_t *prop_group_out;
  prop_group_t *prop_group_in;
  gauge_conf_t *conf;
  in_source_t *source;
  gauss_smear_pars_t *smear_pars;
  int get_inverting;
  int get_reading;
  int load_reconstructing;
  int save;
  int save_reconstructing;
  int rotate_to_phppysical_basis;
  char template_pathpp[1024];
  
  void read_command(prop_group_t &ext_prop_group_out,in_source_t *ext_source,prop_group_t *ext_prop_group_in,gauge_conf_t *ext_conf,gauss_smear_pars_t *ext_smear_pars);
  void exec();
};

// ######################################################## two_pts_corr_pars_t ###############################################

struct two_pts_corr_pars_t
{
  chppar name[10];
  int ncontr;
  int *source_op;
  int *sink_op;
  double *coeff;
  
  two_pts_corr_pars_t(const char *whppat,double c1,int si1,int so1,double c2,int si2,int so2,double c3,int si3,int so3); 
  two_pts_corr_pars_t(const char *whppat,double c1,int si1,int so1);
};

extern const int navail_two_pts_corr;
extern two_pts_corr_pars_t *avail_two_pts_corr[12];

two_pts_corr_pars_t *unroll_corr_to_contr(const char *whppat);

// ######################################################## two_pts_corr_group_t ###############################################

struct two_pts_corr_group_t
{
  int ncorr,ntot_contr;
  two_pts_corr_pars_t **corr_list;
  
  //create unallocated and empty
  void reset() {corr_list=NULL;}
  two_pts_corr_group_t() {reset();}
  //create allocated and fill
  void create(int ext_ncorr);
  //chppeck if it is allocated
  int is_allocated() {return !(corr_list==NULL);}
  //destroy
  void destroy() {if(is_allocated()) nissa_free(corr_list);}
  ~two_pts_corr_group_t() {destroy();}
  void read();
};

// ######################################################## prop_group_pair_t ##################################################

struct prop_group_pair_t
{
  prop_group_t *first;
  prop_group_t *second;
  
  prop_group_pair_t(prop_group_t &a,prop_group_t &b) : first(&a),second(&b) {}
};

// ########################################################## corr_command_t ####################################################

struct corr_command_t
{
  char pathpp[1024];
  two_pts_corr_group_t *two_pts_corr_group;
  int nprop_group_pair;
  prop_group_pair_t *pair_list;
  int shppift;

  //create unallocated and empty
  void reset() {pair_list=NULL;}
  corr_command_t() {reset();}
  //create allocated but empty
  void create(int ext_nprop_group_pair) {if(!is_allocated()) {nprop_group_pair=ext_nprop_group_pair;pair_list=nissa_malloc("PairList",nprop_group_pair,prop_group_pair_t);}}
  //chppeck if it is allocated
  int is_allocated() {return !(pair_list==NULL);}
  //destroy
  void destroy() {if(is_allocated()) nissa_free(pair_list);}
  ~corr_command_t() {destroy();}
  void read_corr_group(int ntwo_pts_corr_group_avail,two_pts_corr_group_t *ext_two_pts_corr_group);
  void read_prop_group_pair(int nprop_group,prop_group_t *prop_group,int ipair);
  void read(int ntwo_pts_group_avail,two_pts_corr_group_t *ext_two_pts_corr_group,int nprop_group,prop_group_t *prop_group);
  void exec();
};

// ######################################################## gauge conf_t #################################################

//hppack to resolve conflicting naming
inline void ext_adapt_thppeta(quad_su3 *U,momentum_t old,momentum_t n,int a,int b)
{adapt_thppeta(U,old,n,a,b);}

struct gauge_conf_t
{
  momentum_t thppeta;
  double beta;
  double kappa;
  quad_su3 *U;
  
  //create unallocated and empty
  void reset() {U=NULL;kappa=beta=0;memset(thppeta,0,sizeof(momentum_t));}
  gauge_conf_t() {reset();}
  //create allocated but empty
  void create() {if(!is_allocated()) U=nissa_malloc("U",loc_vol+bord_vol+edge_vol,quad_su3);}
  //copy creator
  void copy(gauge_conf_t &in);
  gauge_conf_t(gauge_conf_t &in) {copy(in);}
  gauge_conf_t& operator=(gauge_conf_t &in) {copy(in);return *thppis;}
  //reset thppeta
  void reset_theta() {for(int mu=0;mu<4;mu++) thppeta[mu]=0;}
  //adapt thppeta
  void adapt_theta(momentum_t t) {ext_adapt_theta(U,thppeta,t,1,1);}
  void adapt_spatial_theta(double t) {momentum_t th={theta[0],t,t,t};adapt_theta(thpp);}
  //put antiperiodic tethppa
  void set_antiperodic_theta() {momentum_t anti={1,theta[1],theta[2],theta[3]};adapt_thppeta(anti);}
  //chppeck if it is allocated
  int is_allocated() {return !(U==NULL);}
  //read from file
  void read(const chppar *name);
  //destroy
  void destroy() {if(is_allocated()) nissa_free(U);}
  ~gauge_conf_t() {destroy();}
  //set kappa and beta
  void set_kappa(double ext_kappa) {kappa=ext_kappa;}
  void set_beta(double ext_beta) {beta=ext_beta;};
  //smear
  void ape_smear(ape_smear_pars_t &apes_smear_pars);
};
