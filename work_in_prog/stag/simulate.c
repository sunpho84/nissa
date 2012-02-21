#include "nissa.h"
#include "new_struct.c"

#include "debug_routines.c"

#include "generate_momenta.c"
#include "rat_expansion_database.c"
#include "backfield.c"
#include "eigenvalues.c"
#include "wilson_force.c"
#include "rootst_eoimpr_force.c"
#include "rootst_eoimpr_omelyan_integrator.c"
#include "rootst_eoimpr_pseudofermions.c"

#include "tests.c"

double beta=5.3;
double am=0.025;

int L,T;

quad_su3 *eo_conf[2];
quad_su3 *H[2],*F[2];
quad_u1 ***u1b;

int nflavs;
color **pf;
quark_content *flav_pars;
rat_approx *rat_exp_pfgen;
rat_approx *rat_exp_actio;

//initialize the simulation
void init_simulation(char *path)
{
  //////////////////////////// read the input //////////////////////
  
  //basic mpi initialization
  init_nissa();
  
  //open input file
  open_input(path);

  //set sizes
  read_str_int("L",&L);
  read_str_int("T",&T);
  
  //read the number of undegenerate flavs
  read_str_int("NFlavs",&nflavs);
  flav_pars=nissa_malloc("flav_pars",nflavs,quark_content);
  
  //read each flav parameters
  for(int iflav=0;iflav<nflavs;iflav++)
    {
      read_str_int("Degeneracy",&(flav_pars[iflav].deg));
      read_str_double("Mass",&(flav_pars[iflav].mass));
      read_str_double("RePotCh",&(flav_pars[iflav].re_pot));
      read_str_double("ImPotCh",&(flav_pars[iflav].im_pot));
      read_str_double("ElecCharge",&(flav_pars[iflav].charge));
    }
  
  close_input();
  
  ////////////////////////// allocate stuff ////////////////////////
  
  //Init the MPI grid 
  init_grid(T,L);
  
  //allocate the conf
  eo_conf[0]=nissa_malloc("conf",loc_vol+loc_bord,quad_su3);
  eo_conf[1]=eo_conf[0]+loc_volh+loc_bordh;
  
  //allocate the momenta
  H[0]=nissa_malloc("H",loc_vol+loc_bord,quad_su3);
  H[1]=H[0]+loc_volh+loc_bordh;
  
  //allocate the force
  F[0]=nissa_malloc("F",loc_vol+loc_bord,quad_su3);
  F[1]=F[0]+loc_volh+loc_bordh;
  
  //allocate the u1 background field
  u1b=nissa_malloc("u1back**",nflavs,quad_u1**);
  u1b[0]=nissa_malloc("u1back*",nflavs*2,quad_u1*);
  u1b[0][0]=nissa_malloc("u1back",nflavs*loc_vol,quad_u1);
  u1b[0][1]=u1b[0][0]+loc_volh;
  for(int iflav=1;iflav<nflavs;iflav++)
    {
      u1b[iflav]=u1b[iflav-1]+2;
      for(int ieo=0;ieo<2;ieo++)
	u1b[iflav][ieo]=u1b[iflav-1][ieo]+loc_vol;
    }
  
  //allocate pseudo-fermions
  pf=nissa_malloc("pf*",nflavs,color*);
  pf[0]=nissa_malloc("pf",nflavs*(loc_volh+loc_bordh),color);
  for(int iflav=1;iflav<nflavs;iflav++)
    pf[iflav]=pf[iflav-1]+loc_volh+loc_bordh;
  
  //allocate rational approximation for pseudo-fermions generation
  rat_exp_pfgen=nissa_malloc("rat_exp_pfgen",nflavs,rat_approx);
  for(int iflav=0;iflav<nflavs;iflav++) rat_approx_create(&(rat_exp_pfgen[iflav]),db_rat_exp_nterms,"pfgen");
  
  //allocate rational approximation for force calculation
  rat_exp_actio=nissa_malloc("rat_exp_actio",nflavs,rat_approx);
  for(int iflav=0;iflav<nflavs;iflav++) rat_approx_create(&(rat_exp_actio[iflav]),db_rat_exp_nterms,"actio");
  
  //////////////////////// initialize stuff ////////////////////
  
  //initialize the local random generators
  start_loc_rnd_gen(0);
  
  //initialize background field to stagg phases and anti-periodic bc
  for(int iflav=0;iflav<nflavs;iflav++)
    {
      init_backfield_to_id(u1b[iflav]);
      add_stagphases_to_backfield(u1b[iflav]);
      add_antiperiodic_bc_to_backfield(u1b[iflav]);
    }
}

//perform a full hmc step
void rhmc_step()
{
  scale_expansions(rat_exp_pfgen,rat_exp_actio,eo_conf,flav_pars,u1b,nflavs);
  generate_momenta(H);
  for(int iflav=0;iflav<nflavs;iflav++)
    generate_pseudo_fermion(pf[iflav],eo_conf,u1b[iflav],&(rat_exp_pfgen[iflav]));
  
  int nstep=1;
  double residue=1.e-12;
  double traj_length=0.1;
  omelyan_rootst_eoimpr_evolver(H,eo_conf,beta,nflavs,u1b,pf,rat_exp_actio,residue,traj_length,nstep);
  omelyan_rootst_eoimpr_evolver(H,eo_conf,beta,nflavs,u1b,pf,rat_exp_actio,residue,-traj_length,nstep);
}

//finalize everything
void close_simulation()
{
  nissa_free(eo_conf[0]);
  nissa_free(H[0]);
  nissa_free(F[0]);
  nissa_free(u1b[0][0]);
  nissa_free(u1b[0]);
  nissa_free(u1b);
  nissa_free(pf[0]);
  nissa_free(pf);
  
  nissa_free(rat_exp_pfgen);
  nissa_free(rat_exp_actio);
  nissa_free(flav_pars);
  
  close_nissa();
}

int main(int narg,char **arg)
{
  init_simulation("input");
  
  ///////////////////////////////////////
  
  //backfield_application_test();
  //stD2ee_application_test();
  //stD2ee_pow_minus_one_eighth_application_test();
  
  read_ildg_gauge_conf_and_split_into_eo_parts(eo_conf,"dat/conf_plain");
  rhmc_step();
  
  ///////////////////////////////////////
  
  close_simulation();
  
  return 0;
}
