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

double beta=5.3;
double am=0.025;

int L,T;

quad_su3 *new_conf[2];
quad_su3 *conf[2];
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
  conf[0]=nissa_malloc("conf",loc_vol+loc_bord,quad_su3);
  conf[1]=conf[0]+loc_volh+loc_bordh;
  new_conf[0]=nissa_malloc("new_conf",loc_vol+loc_bord,quad_su3);
  new_conf[1]=new_conf[0]+loc_volh+loc_bordh;
  
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
  
  //initialize background field to id
  for(int iflav=0;iflav<nflavs;iflav++)
    {
      init_backfield_to_id(u1b[iflav]);
    }
}

//copy an e/o split conf
void eo_conf_copy(quad_su3 **dest,quad_su3 **source)
{
  for(int eo=0;eo<2;eo++)
    memcpy(dest[eo],source[eo],loc_volh*sizeof(quad_su3));
}

//perform a full hmc step
void rhmc_step(quad_su3 **out_conf,quad_su3 **in_conf)
{
  //copy the old conf into the new
  eo_conf_copy(out_conf,in_conf);
  
  //add the phases
  addrem_stagphases_to_eo_conf(out_conf);
  
  //generate the appropriate expansion of rational approximations
  scale_expansions(rat_exp_pfgen,rat_exp_actio,conf,flav_pars,u1b,nflavs);
  
  //create the momenta
  generate_momenta(H);
  
  //create pseudo-fermions
  for(int iflav=0;iflav<nflavs;iflav++)
    generate_pseudo_fermion(pf[iflav],out_conf,u1b[iflav],&(rat_exp_pfgen[iflav]));
  
  int nstep=13;
  double residue=1.e-12;
  double traj_length=1;
  omelyan_rootst_eoimpr_evolver(H,out_conf,beta,nflavs,u1b,pf,rat_exp_actio,residue,traj_length,nstep);
  omelyan_rootst_eoimpr_evolver(H,out_conf,beta,nflavs,u1b,pf,rat_exp_actio,residue,-traj_length,nstep);
  
  //remove the phases
  addrem_stagphases_to_eo_conf(out_conf);
}

//finalize everything
void close_simulation()
{
  nissa_free(new_conf[0]);
  nissa_free(conf[0]);
  
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

//compare an eo_conf with a saved one
void check_eo_conf(quad_su3 **eo_conf,char *path)
{
  //allocate
  quad_su3 *temp[2];
  temp[0]=nissa_malloc("temp",loc_vol,quad_su3);
  temp[1]=temp[0]+loc_volh;
  
  //read
  master_printf("Debug, reading conf after updating\n");
  read_ildg_gauge_conf_and_split_into_eo_parts(temp,path);
  
  //compute the norm
  double n2=0;
  for(int eo=0;eo<2;eo++)
    for(int ivol=0;ivol<loc_volh;ivol++)
      for(int mu=0;mu<4;mu++)
	for(int ic1=0;ic1<3;ic1++)
          for(int ic2=0;ic2<3;ic2++)
            for(int ri=0;ri<2;ri++)
              {
		double a=temp[eo][ivol][mu][ic1][ic2][ri]-eo_conf[eo][ivol][mu][ic1][ic2][ri];
		n2+=a*a;
              }
  n2/=loc_vol*4*9;
  n2=sqrt(n2);
  
  //check
  master_printf("Total conf norm diff: %lg\n",n2);
  if(n2>1.e-7) crash("conf updating failed");
  
  //delocate
  nissa_free(temp[0]);
}

int main(int narg,char **arg)
{
  init_simulation("input");
  
  ///////////////////////////////////////
  
  read_ildg_gauge_conf_and_split_into_eo_parts(conf,"dat/conf_plain");
  rhmc_step(new_conf,conf);
  
  //debug
  //check_eo_conf(new_conf,"dat/final_conf");
  check_eo_conf(new_conf,"dat/conf_plain");
  
  ///////////////////////////////////////
  
  close_simulation();
  
  return 0;
}
