#pragma once

#include "momenta_action.c"

//compute quark action for a set of quark
double rootst_eoimpr_quark_action(quad_su3 **eo_conf,int nfl,quad_u1 ***u1b,color **pf,rat_approx *appr,double residue)
{  
  //allocate chi
  color *chi_e=nissa_malloc("chi_e",loc_volh+loc_bordh,color);
  
  //quark action
  double loc_action=0;
  for(int ifl=0;ifl<nfl;ifl++)
    {
      //compute chi with background field
      add_backfield_to_conf(eo_conf,u1b[ifl]);
      summ_src_and_all_inv_stD2ee_cgmm2s(chi_e,pf[ifl],eo_conf,appr,1000000,residue,residue,0);
      rem_backfield_from_conf(eo_conf,u1b[ifl]);
      
      //compute scalar product
      nissa_loc_volh_loop(ivol)
	for(int ic=0;ic<3;ic++)
	  loc_action+=chi_e[ivol][ic][0]*pf[ifl][ivol][ic][0]+chi_e[ivol][ic][1]*pf[ifl][ivol][ic][1];
    }
  
  //global reducton
  double glb_action;
  MPI_Allreduce(&loc_action,&glb_action,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  
  //free
  nissa_free(chi_e);
  
  return glb_action;
}

//Compute the total action of the rooted staggered e/o improved theory.
//Passed conf must NOT contain the backfield.
double full_rootst_eoimpr_action(quad_su3 **eo_conf,double beta,quad_su3 **H,int nfl,quad_u1 ***u1b,color **pf,rat_approx *appr,double residue)
{
  master_printf("Computing action\n");

  //compute the three parts of the action
  double quark_action=rootst_eoimpr_quark_action(eo_conf,nfl,u1b,pf,appr,residue);
  
  //gauge action
  double gluon_action=beta*6*(1+global_plaquette_eo_conf(eo_conf[0],eo_conf[1]))*glb_vol;
  
  //momenta action
  double mom_action=momenta_action(H);
  
  //*debug
  master_printf("Q%.18lg\n",quark_action);
  master_printf("G%.18lg\n",gluon_action);
  master_printf("M%.18lg\n",mom_action);
  //*/
  
  return quark_action+gluon_action+mom_action;
}
