#include "nissa.h"
#include "new_struct.c"

#include "debug_routines.c"

#include "rat_expansion_database.c"
#include "backfield.c"
#include "eigenvalues.c"
#include "rootst_eo_impr_omelyan_integrator.c"
#include "rootst_eo_impr_force.c"
#include "wilson_force.c"

#include "tests.c"

double beta=5.38;
double am=0.025;

int L,T;
quad_su3 *eo_conf[2];
color *pf[2];
int npf=2;
quad_su3 *H[2],*F[2];
quad_u1 ***u1b;
rat_approx *rat_exp_pfgen;

int nquarks;
quark_content *quark_pars;

//initialize the simulation
void init_simulation(char *path)
{
  //basic mpi initialization
  init_nissa();
  
  //open input file
  open_input(path);

  //set sizes
  read_str_int("L",&L);
  read_str_int("T",&T);
  
  //Init the MPI grid 
  init_grid(T,L);
  
  //initialize the local random generators
  start_loc_rnd_gen(0);
  
  //allocate the conf
  eo_conf[0]=nissa_malloc("conf",loc_vol+loc_bord,quad_su3);
  eo_conf[1]=eo_conf[0]+(loc_vol+loc_bord)/2;
  
  //allocate the momenta
  H[0]=nissa_malloc("H",loc_vol+loc_bord,quad_su3);
  H[1]=H[0]+(loc_vol+loc_bord)/2;
  
  //allocate the force
  F[0]=nissa_malloc("F",loc_vol+loc_bord,quad_su3);
  F[1]=F[0]+(loc_vol+loc_bord)/2;
  
  //allocate the u1 background field
  u1b=nissa_malloc("u1back**",nquarks,quad_u1**);
  u1b[0]=nissa_malloc("u1back*",nquarks*2,quad_u1*);
  u1b[0][0]=nissa_malloc("u1back",nquarks*loc_vol,quad_u1);
  u1b[0][1]=u1b[0][0]+loc_vol/2;
  for(int iquark=1;iquark<nquarks;iquark++)
    {
      u1b[iquark]=u1b[iquark-1]+2;
      for(int ieo=0;ieo<2;ieo++)
	u1b[iquark][ieo]=u1b[iquark-1][ieo]+loc_vol;
    }
  
  //allocate pseudo-fermions
  for(int iq=0;iq<2;iq++)
    pf[iq]=nissa_malloc("pf",(loc_vol+loc_bord)/2,color);
  
  //read the conf
  read_ildg_conf_and_split_into_eo_parts(eo_conf,"dat/conf_plain");
  
  //initialize background field to stagg phases and anti-periodic bc
  for(int iquark=0;iquark<nquarks;iquark++)
    {
      init_backfield_to_id(u1b[iquark]);
      add_stagphases_to_backfield(u1b[iquark]);
      add_antiperiodic_bc_to_backfield(u1b[iquark]);
    }
  
  //read the number of undegenerate quarks
  read_str_int("nquarks",&nquarks);
  quark_pars=nissa_malloc("quark_pars",nquarks,quark_content);
  
  //read each quark parameters
  for(int iquark=0;iquark<nquarks;iquark++)
    {
      read_str_int("Degeneracy",&(quark_pars[iquark].deg));
      read_str_double("Mass",&(quark_pars[iquark].mass));
      read_str_double("RePotCh",&(quark_pars[iquark].re_pot));
      read_str_double("ImPotCh",&(quark_pars[iquark].im_pot));
      read_str_double("ElecCharge",&(quark_pars[iquark].charge));
    }
  
  //allocate rational approximation for pseudo-fermions generation
  rat_exp_pfgen=nissa_malloc("rat_exp_pfgen",nquarks,rat_approx);
  for(int iquark=0;iquark<nquarks;iquark++) rat_approx_create(&(rat_exp_pfgen[iquark]),db_rat_exp_nterms,"pfgen");
  
  close_input();
}

//generate momenta using guassian hermitean matrix generator
void generate_momenta()
{
  for(int eo=0;eo<2;eo++)
    for(int ivol=0;ivol<loc_vol/2;ivol++)
      for(int mu=0;mu<4;mu++)
	herm_put_to_gauss(H[eo][ivol][mu],&(loc_rnd_gen[ivol]),1);
}

//perform a full hmc step
void rhmc_step()
{
  scale_expansions(rat_exp_pfgen,eo_conf,quark_pars,u1b,nquarks);
  generate_momenta();
  generate_pseudo_fermions();
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
  
  for(int iq=0;iq<2;iq++)
    nissa_free(pf[iq]);
  
  nissa_free(rat_exp_pfgen);
  nissa_free(quark_pars);
  
  close_nissa();
}

int main(int narg,char **arg)
{
  init_simulation("input");
  
  ///////////////////////////////////////
  
  read_ildg_conf_and_split_into_eo_parts(eo_conf_comp,"dat/conf_plain");
    
  ///////////////////////////////////////
  
  close_simulation();
  
  return 0;
}
