#include "../new_types/new_types_definitions.h"
#include "../new_types/complex.h"
#include "../base/global_variables.h"
#include "../base/macros.h"
#include "../base/random.h"
#include "../base/routines.h"
#include "../base/vectors.h"
#include "../geometry/geometry_eo.h"
#include "../hmc/backfield.h"
#include "../inverters/staggered/cg_invert_stD.h"

//get a propagator
void get_propagator(color **prop,quad_su3 **conf,quad_u1 **u1b,double m,double residue,color **source)
{
  addrem_stagphases_to_eo_conf(conf);
  add_backfield_to_conf(conf,u1b);
  inv_stD_cg(prop,conf,m,10000,5,residue,source);
  rem_backfield_from_conf(conf,u1b);
  addrem_stagphases_to_eo_conf(conf);
}

//compute the chiral condensate
void chiral_condensate(complex cond,quad_su3 **conf,quad_u1 **u1b,double m,double residue)
{
  //allocate
  color *rnd[2]={nissa_malloc("rnd_EVN",loc_volh+bord_volh,color),nissa_malloc("rnd_ODD",loc_volh+bord_volh,color)};
  color *chi[2]={nissa_malloc("chi_EVN",loc_volh+bord_volh,color),nissa_malloc("chi_ODD",loc_volh+bord_volh,color)};
  
  //generate the source and the propagator
  generate_fully_undiluted_source(rnd,RND_Z4,-1);
  get_propagator(rnd,conf,u1b,m,residue,rnd);
  
  //compute the condensate
  complex loc_cond={0,0};
  for(int par=0;par<2;par++)
    nissa_loc_volh_loop(ivol)
      for(int icol=0;icol<3;icol++)
	complex_summ_the_conj2_prod(loc_cond,chi[par][ivol][icol],rnd[par][ivol][icol]);
  
  //global reduction
  glb_reduce_complex(cond,loc_cond);
  
  //add normalization: 1/4vol
  complex_prodassign_double(cond,1.0/(4*glb_vol));
  
  //free
  for(int par=0;par<2;par++)
    {
      nissa_free(rnd[par]);
      nissa_free(chi[par]);
    }
}

//measure chiral cond
void measure_chiral_cond(quad_su3 **conf,theory_pars_type &theory_pars,int iconf)
{
  FILE *file=open_file(theory_pars.chiral_cond_pars.path,(iconf==0)?"w":"a");

  master_fprintf(file,"%d",iconf);
  
  //measure the condensate for each quark
  for(int iflav=0;iflav<theory_pars.nflavs;iflav++)
    {
      complex cond={0,0};
      
      //loop over hits
      int nhits=theory_pars.chiral_cond_pars.nhits;
      for(int hit=0;hit<nhits;hit++)
        {
          verbosity_lv1_master_printf("Evaluating chiral condensate for flavor %d/%d, nhits %d/%d\n",iflav+1,theory_pars.nflavs,hit+1,nhits);
          
          //compute and summ
          complex temp;
          chiral_condensate(temp,conf,theory_pars.backfield[iflav],theory_pars.quark_content[iflav].mass,theory_pars.chiral_cond_pars.residue);
          complex_summ_the_prod_double(cond,temp,1.0/nhits);
        }
      
      master_fprintf(file,"\t%16.16lg",cond[0]);
    }

  master_fprintf(file,"\n");
  
  if(rank==0) fclose(file);
}
