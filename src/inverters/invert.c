#pragma once

//follow inclusion below
#include "cgmms_invert_tmQ2_common.c"
#include "cgmms_invert_stD2ee_common.c"

#ifdef BGP

#include "cg_invert_tmQ2_bgp.c"
#include "cgmms_invert_tmQ2_bgp.c"

#include "cgmms_invert_stD2ee_bgp.c"

#else

#include "cg_invert_tmDeoimpr_portable.c"
#include "cg_invert_tmQ2_portable.c"
#include "cgmms_invert_tmQ2_portable.c"

#include "cgmms_invert_stD2ee_portable.c"

#endif

#include "cg_invert_tmQ2_common.c"

////////////////////////////////////// TWISTED MASS LIGHT DEGENERATE INVERTERS ///////////////////////////////

void inv_tmQ2_cgmms(spincolor **sol,quad_su3 *conf,double kappa,double *m,int nmass,int niter,double st_res,double st_minres,int st_crit,spincolor *source)
{inv_tmQ2_cgmms_RL(sol,conf,kappa,m,nmass,niter,st_res,st_minres,st_crit,0,source);}
void inv_tmQ2_cgmms_left(spincolor **sol,quad_su3 *conf,double kappa,double *m,int nmass,int niter,double st_res,double st_minres,int st_crit,spincolor *source)
{inv_tmQ2_cgmms_RL(sol,conf,kappa,m,nmass,niter,st_res,st_minres,st_crit,1,source);}

void inv_tmDQ_cgmms_RL(spincolor **sol,quad_su3 *conf,double kappa,double *m,int nmass,int niter,double st_res,double st_minres,int st_crit,int RL,spincolor *source)
{
  //put the g5
  nissa_loc_vol_loop(ivol) for(int id1=2;id1<4;id1++) for(int ic1=0;ic1<3;ic1++) for(int ri=0;ri<2;ri++) source[ivol][id1][ic1][ri]*=-1;
  set_borders_invalid(source);
  inv_tmQ2_cgmms_RL(sol,conf,kappa,m,nmass,niter,st_res,st_minres,st_crit,RL,source);
  nissa_loc_vol_loop(ivol) for(int id1=2;id1<4;id1++) for(int ic1=0;ic1<3;ic1++) for(int ri=0;ri<2;ri++) source[ivol][id1][ic1][ri]*=-1;
  set_borders_invalid(source);
}

void inv_tmDQ_cgmms(spincolor **sol,quad_su3 *conf,double kappa,double *m,int nmass,int niter,double st_res,double st_minres,int st_crit,spincolor *source)
{inv_tmDQ_cgmms_RL(sol,conf,kappa,m,nmass,niter,st_res,st_minres,st_crit,0,source);}
void inv_tmDQ_cgmms_left(spincolor **sol,quad_su3 *conf,double kappa,double *m,int nmass,int niter,double st_res,double st_minres,int st_crit,spincolor *source)
{inv_tmDQ_cgmms_RL(sol,conf,kappa,m,nmass,niter,st_res,st_minres,st_crit,1,source);}

////////////////////////////////////////////////// full frontend //////////////////////////////////////////////

//invert a set of propagators using the passed source
//the output is stored in twisted basis, assuming that prop=su3spinspin[2][nmass][>=loc_vol]
void compute_su3spinspin_propagators_multi_mass(su3spinspin ***prop,quad_su3 *conf,double kappa,double *mass,int nmass,int niter_max,double stopping_residue,double minimal_residue,int stopping_criterion,su3spinspin *source)
{
  //allocate temporary source
  spincolor *temp_source=nissa_malloc("temp_source",loc_vol+loc_bord,spincolor);
  //allocate temp_vec
  spincolor *temp_vec[2];
  temp_vec[0]=nissa_malloc("temp_vec[0]",loc_vol+loc_bord,spincolor);
  temp_vec[1]=nissa_malloc("temp_vec[1]",loc_vol+loc_bord,spincolor);
  //allocate nmass spincolors, for the cgmms solutions
  spincolor **cgmms_solution;
  cgmms_solution=nissa_malloc("cgmms_solution",nmass,spincolor*);
  for(int imass=0;imass<nmass;imass++) cgmms_solution[imass]=nissa_malloc("cgmms_solution[imass]",loc_vol+loc_bord,spincolor);
  
  //loop over the source dirac and color index
  for(int id=0;id<4;id++)
    for(int ic=0;ic<3;ic++)
      {
        nissa_loc_vol_loop(ivol)
          get_spincolor_from_su3spinspin(temp_source[ivol],source[ivol],id,ic);
        set_borders_invalid(temp_source);
	
        double init_time=take_time();
        inv_tmDQ_cgmms(cgmms_solution,conf,kappa,mass,nmass,niter_max,stopping_residue,minimal_residue,stopping_criterion,temp_source);
        master_printf("Finished the inversion of Q2, dirac index %d, color %d in %g sec\n",id,ic,take_time()-init_time);
        
        //reconstruct the doublet
        for(int imass=0;imass<nmass;imass++)
          {
            reconstruct_tm_doublet(temp_vec[0],temp_vec[1],conf,kappa,mass[imass],cgmms_solution[imass]);
            
            //convert the id-th spincolor into the su3spinspin
            for(int r=0;r<2;r++)
	      {
		nissa_loc_vol_loop(ivol)       
		  put_spincolor_into_su3spinspin(prop[r][imass][ivol],temp_vec[r][ivol],id,ic);
		set_borders_invalid(prop[r]);
	      }
            master_printf("Mass %d (%g) reconstructed \n",imass,mass[imass]);
          }
      }
  
  //free temp vec
  nissa_free(temp_vec[1]);nissa_free(temp_vec[0]);
  nissa_free(temp_source);
  for(int imass=0;imass<nmass;imass++) nissa_free(cgmms_solution[imass]);
  nissa_free(cgmms_solution);
}
