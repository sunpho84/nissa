#ifndef _PROP_HPP
#define _PROP_HPP

//invert on top of a source, putting all needed for the appropriate quark
inline void get_qprop(PROP_TYPE *out,PROP_TYPE *in,int imass,bool r)
{
  //these are the ways in which Dirac operator rotates - propagator is opposite, see below
#ifdef POINT_SOURCE_VERSION
  for(int ic=0;ic<NCOL;ic++)
#endif
    for(int id=0;id<4;id++)
      { 
	//read the source out
#ifdef POINT_SOURCE_VERSION
	get_spincolor_from_su3spinspin(temp_source,in,id,ic);
#else
	get_spincolor_from_colorspinspin(temp_source,in,id);
#endif
	
	//rotate the source index - the propagator rotate AS the sign of mass term
	if(!pure_wilson) safe_dirac_prod_spincolor(temp_source,(tau3[r]==-1)?&Pminus:&Pplus,temp_source);
	
	//invert
	inv_time-=take_time();
	if(!pure_wilson) inv_tmD_cg_eoprec_eos(temp_solution,NULL,conf,kappa,tau3[r]*qmass[imass],100000,residue[imass],temp_source);
	else             inv_tmD_cg_eoprec_eos(temp_solution,NULL,conf,qkappa[imass],0,100000,residue[imass],temp_source);
	ninv_tot++;inv_time+=take_time();
	
	//rotate the sink index
	if(!pure_wilson) safe_dirac_prod_spincolor(temp_solution,(tau3[r]==-1)?&Pminus:&Pplus,temp_solution);
	
	//put the output on place
#ifdef POINT_SOURCE_VERSION
	master_printf("  finished the inversion dirac index %d, color %d\n",id,ic);
	put_spincolor_into_su3spinspin(out,temp_solution,id,ic);
#else
	master_printf("  finished the inversion dirac index %d\n",id);
	put_spincolor_into_colorspinspin(out,temp_solution,id);
#endif
      }
}

#endif
