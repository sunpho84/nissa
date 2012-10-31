#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include "../../src/nissa.h"

struct mass_res_group_t;
struct theta_group_t;
struct ape_smear_pars_t;
struct gauss_smear_pars_t;
struct source_t;
struct prop_group_t;
struct prop_group_command_t;
struct gauge_conf_t;

extern char base_out_folder[1024];

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
  //check if it is allocated
  int is_allocated() {return !(mass==residues && mass==NULL);}
  //destroy
  void destroy() {if(is_allocated()){nissa_free(mass);nissa_free(residues);}}
  ~mass_res_group_t() {destroy();}
  //read from input file
  void read() {destroy();int n;read_str_int("NMassRes",&n);create(n);for(int i=0;i<n;i++){read_double(&mass[i]);read_double(&residues[i]);}}
  //print to screen
  void printf() {master_printf("List of masses and residues\n");for(int i=0;i<nmass;i++){master_printf("%lg %lg\n",mass[i],residues[i]);}}
};

// ######################################################### theta_group_t ###################################################

struct theta_group_t
{
  int ntheta;
  double *theta;
  
  //create unallocated and empty
  theta_group_t() {theta=NULL;}
  //create allocated but empty
  void create(int n) {destroy();ntheta=n;theta=nissa_malloc("theta",n,double);}
  theta_group_t(int n) {create(n);}
  //check if it is allocated
  int is_allocated() {return !(theta==NULL);}
  //destroy
  void destroy() {if(is_allocated()) nissa_free(theta);}
  ~theta_group_t() {destroy();}
  //read from input file
  void read() {destroy();int n;read_str_int("NTheta",&n);create(n);for(int i=0;i<n;i++) read_double(&theta[i]);}
  //print to screen
  void printf() {master_printf("List of thetas\n");for(int i=0;i<ntheta;i++) master_printf("%lg\n",theta[i]);}
};

// ###################################################### ape_smear_pars_t ###############################################

struct ape_smear_pars_t
{
  double alpha;
  int niter;
  
  //create with default values
  ape_smear_pars_t(int ext_niter=20,double ext_alpha=0.5){niter=ext_niter;alpha=ext_alpha;}
  //read from input file
  void read() {read_str_double("ApeAlpha",&alpha);read_str_int("ApeNiter",&niter);}
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
  //check if it is allocated
  int is_allocated() {return !(coeff==NULL);}
  //destroy
  void destroy() {if(is_allocated()) {nissa_free(coeff);nissa_free(expnt);}}
  ~gauss_smear_pars_t() {destroy();}
  //read from input file
  void read();
  //set kappa manually
  void set_kappa(double ext_kappa) {kappa=ext_kappa;}
};

// ########################################################## source_t ###################################################

struct source_t
{
  colorspinspin *eta;
  
  //create unallocated and empty
  void reset() {eta=NULL;}
  source_t() {reset();}
  //create allocated but empty
  void create() {if(!is_allocated()) eta=nissa_malloc("eta",loc_vol+bord_vol,colorspinspin);}
  //copy and assignement
  void copy(const source_t &in) {create();vector_copy(eta,in.eta);}
  source_t(const source_t &in) {copy(in);}
  source_t& operator=(const source_t &in) {copy(in);return *this;}
  //check if it is allocated
  int is_allocated() {return !(eta==NULL);}
  //destroy
  void destroy() {if(is_allocated()) nissa_free(eta);}
  ~source_t() {destroy();}
  //read from file
  void read(const char *name) {create();read_colorspinspin(eta,name,NULL);}
  //write to file
  void write(const char *name) {write_colorspinspin(name,eta,64);}
  //fill a particular timeslice
  void fill(rnd_type noise_type,int twall) {create();generate_spindiluted_source(eta,noise_type,twall);}
  //smear using a gaussian smearing operator
  void smear(gauge_conf_t &conf,gauss_smear_pars_t &pars);
};

// ####################################################### prop_group_t ###################################################

enum TMR{R_ZERO,R_ONE,R_BOTH};

struct prop_group_t
{
  mass_res_group_t *mass_res;
  theta_group_t *theta;
  TMR which_r;
  colorspinspin **S;
  
  //create unallocated and empty
  void reset(){S=NULL;}
  prop_group_t() {reset();}
  //create allocated but empty
  void create(theta_group_t &t,mass_res_group_t &m,TMR r);
  prop_group_t(theta_group_t &t,mass_res_group_t &m,TMR r) {create(t,m,r);}
  //check if it is allocated
  int is_allocated() {return !(S==NULL);}
  //check the validity of the passed args
  void check_itheta_mass_r(int itheta,int imass,int r);
  //return number of theta,masses and r
  void get_ntheta_mass_r(int &ntheta,int &nmass,int &nr);
  //return the id of the passed args
  int iprop(int itheta,int imass,int r);
  //compute the number of propagators
  int nprop();
  //destroy
  void destroy() {if(is_allocated()) {for(int i=0,n=nprop();i<n;i++) nissa_free(S[i]);nissa_free(S);}}
  ~prop_group_t() {destroy();}
  //invert
  void get_inverting(source_t &source,gauge_conf_t &gauge_conf,int rotate_to_physical_basis);
  //smear
  void get_smearing(gauge_conf_t &conf,gauss_smear_pars_t &pars,prop_group_t &in);
  void smear(gauge_conf_t &conf,gauss_smear_pars_t &pars);
  //read the parameters
  void read_pars(int ntheta_group,theta_group_t *t,int nmass_res_group,mass_res_group_t *m);
  //read the group
  void get_reading(const char *template_path,gauge_conf_t &conf,int load_reconstructing,int load_rotating_to_physical_basis);
  //write the group
  void write(const char *template_path,int save_reconstructing,int is_rotated,gauge_conf_t &conf);
};

// ################################################### prop_group_command_t ###############################################

struct prop_group_command_t
{
  prop_group_t *prop_group_out;
  prop_group_t *prop_group_in;
  gauge_conf_t *conf;
  source_t *source;
  gauss_smear_pars_t *smear_pars;
  int get_inverting;
  int get_reading;
  int load_reconstructing;
  int save;
  int save_reconstructing;
  int rotate_to_physical_basis;
  char template_path[1024];
  
  void read_command(prop_group_t &ext_prop_group_out,source_t *ext_source,prop_group_t *ext_prop_group_in,gauge_conf_t *ext_conf,gauss_smear_pars_t *ext_smear_pars);
  void exec();
};

// ######################################################## two_pts_contr_pars_t ###############################################

struct two_pts_contr_pars_t
{
  int source_op;
  int sink_op;
  double coeff;
  int starting;
  
  two_pts_contr_pars_t(int ext_starting,int ext_sink_op,int ext_source_op,double ext_coeff) : starting(ext_starting),sink_op(ext_sink_op),source_op(ext_source_op),coeff(ext_coeff) {}
};

// ######################################################## two_pts_corr_group_t ###############################################

struct two_pts_corr_group_t
{
  std::vector<two_pts_contr_pars_t> contr_list;
  std::vector<std::string> corr_name;
  
  void add_corr(const char *what);
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
  char path[1024];
  two_pts_corr_group_t *two_pts_corr_group;
  std::vector<prop_group_pair_t> pair_list;
  int shift;

  void read_corr_group(int ntwo_pts_corr_group_avail,two_pts_corr_group_t *ext_two_pts_corr_group);
  void read_prop_group_pair(int nprop_group,prop_group_t *prop_group);
  void read(int ntwo_pts_group_avail,two_pts_corr_group_t *ext_two_pts_corr_group,int nprop_group,prop_group_t *prop_group);
  void exec();
};

// ######################################################## gauge conf_t #################################################

//hack to resolve conflicting naming
inline void ext_adapt_theta(quad_su3 *U,momentum_t old,momentum_t n,int a,int b)
{adapt_theta(U,old,n,a,b);}

struct gauge_conf_t
{
  momentum_t theta;
  double beta;
  double kappa;
  quad_su3 *U;
  
  //create unallocated and empty
  void reset() {U=NULL;kappa=beta=0;memset(theta,0,sizeof(momentum_t));}
  gauge_conf_t() {reset();}
  //create allocated but empty
  void create() {if(!is_allocated()) U=nissa_malloc("U",loc_vol+bord_vol+edge_vol,quad_su3);}
  //copy creator
  void copy(gauge_conf_t &in);
  gauge_conf_t(gauge_conf_t &in) {copy(in);}
  gauge_conf_t& operator=(gauge_conf_t &in) {copy(in);return *this;}
  //reset theta
  void reset_theta() {for(int mu=0;mu<4;mu++) theta[mu]=0;}
  //adapt theta
  void adapt_theta(momentum_t t) {ext_adapt_theta(U,theta,t,1,1);}
  void adapt_spatial_theta(double t) {momentum_t th={theta[0],t,t,t};adapt_theta(th);}
  //put antiperiodic tetha
  void set_antiperodic_theta() {momentum_t anti={1,theta[1],theta[2],theta[3]};adapt_theta(anti);}
  //check if it is allocated
  int is_allocated() {return !(U==NULL);}
  //read from file
  void read(const char *name);
  //destroy
  void destroy() {if(is_allocated()) nissa_free(U);}
  ~gauge_conf_t() {destroy();}
  //set kappa and beta
  void set_kappa(double ext_kappa) {kappa=ext_kappa;}
  void set_beta(double ext_beta) {beta=ext_beta;};
  //smear
  void ape_smear(ape_smear_pars_t &apes_smear_pars);
};
