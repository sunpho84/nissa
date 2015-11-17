#ifndef _PROP_HPP
#define _PROP_HPP

#include "nissa.hpp"

#include "conf.hpp"
#include "pars.hpp"

namespace nissa
{
  int ninv_tot=0;
  double inv_time=0;
  
  inline void get_qprop(spincolor *out,spincolor *in,int imass,bool r)
  {
    //rotate the source index - the propagator rotate AS the sign of mass term
    if(!pure_wilson) safe_dirac_prod_spincolor(in,(tau3[r]==-1)?&Pminus:&Pplus,in);
    
    //invert
    inv_time-=take_time();
    if(!pure_wilson) inv_tmD_cg_eoprec_eos(out,NULL,conf,kappa,tau3[r]*qmass[imass],100000,residue[imass],in);
    else             inv_tmD_cg_eoprec_eos(out,NULL,conf,qkappa[imass],0,100000,residue[imass],in);
    ninv_tot++;inv_time+=take_time();
    
    //rotate the sink index
    if(!pure_wilson) safe_dirac_prod_spincolor(out,(tau3[r]==-1)?&Pminus:&Pplus,out);
  }
  
  //invert on top of a source, putting all needed for the appropriate quark
  inline void get_qprop(PROP_TYPE *out,PROP_TYPE *in,int imass,bool r)
  {
    spincolor *temp_source;
    spincolor *temp_solution;
    
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

	  get_qprop(temp_solution,temp_source,imass,r);
	  
	  //put the output on place
#ifdef POINT_SOURCE_VERSION
	  master_printf("  finished the inversion dirac index %d, color %d\n",id,ic);
	  put_spincolor_into_su3spinspin(out,temp_solution,id,ic);
#else
	  master_printf("  finished the inversion dirac index %d\n",id);
	  put_spincolor_into_colorspinspin(out,temp_solution,id);
#endif
	}
    
    nissa_free(temp_source);
    nissa_free(temp_solution);
  }
}

#endif
