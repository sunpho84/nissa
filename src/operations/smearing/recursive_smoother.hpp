#ifndef _RECURSIVE_SMOOTHER_HPP
#define _RECURSIVE_SMOOTHER_HPP

#include "base/vectors.hpp"
#include "geometry/geometry_lx.hpp"
#include "new_types/su3.hpp"
#include "operations/smearing/smooth.hpp"

#include <vector>

namespace nissa
{
  //holds the data to recursive smooth, generalizing  from 1302.5246
  struct recursive_smoother_t
  {
    struct smooth_lev_t
    {
      quad_su3 *conf; //stored configuration
      int nsmooth; //current smooth in terms of the fundamental
      int epsb; //length of the smoothing in terms of the fundamental, for this level
    };
    
    int mb; //block size
    int nb; //number of blocks at each level
    smooth_pars_t &smoother; //smoother embedded
    
    //level of smoothed conf
    std::vector<smooth_lev_t> levs;
    
    //number of levels
    int nl()
    {return levs.size();}
    
    recursive_smoother_t(quad_su3 **ori_conf,int ns,smooth_pars_t &smoother) : smoother(smoother)
    {
      //take the number of smooth levels
      int m=smoother.nsmooth();
      
      //check
      if(ns<0 or ns>m) crash("ns=%d must be in the range [%d,%d]",0,m);
      
      //set the optimal number of blocks
      nb=pow(m,1.0/(ns+1.0));
      
      //resize
      levs.resize(ns+1);
      for(int i=0;i<=ns;i++)
	{
	  //allocate
	  levs[i].conf=nissa_malloc("conf",loc_vol+bord_vol,quad_su3);
	  //take the first
	  vector_copy(levs[i].conf,ori_conf);
	  //set smooth to 0
	  levs[i].nsmooth=0;
	}
      
      //set units of smooth, recursively
      levs.back().epsb=1;
      for(int is=ns-1;is>=0;is++)
	levs[is].epsb=levs[is+1].epsb*mb;
    }
    
    //find the level which has the closest smaller nsmooth
    int find_closest_smaller_nsmooth(int nsmooth)
    {
      master_printf(" searching the level closest to %d\n",nsmooth);
      
      int iclosest_smooth=0;
      for(int is=1;is<nl();is++)
	{
	  int nclosest_smooth=levs[iclosest_smooth].nsmooth;
	  int lev_ncur_smooth=levs[is].nsmooth;
	  master_printf(" level %d: %d, current closest: %d, ncloses_smooth: %d\n",is,lev_ncur_smooth,iclosest_smooth,nclosest_smooth);
	  if(lev_ncur_smooth<nsmooth and lev_ncur_smooth>nclosest_smooth) iclosest_smooth=is;
	}
      
      return iclosest_smooth;
    }
    
    //copy level
    void copy_level(int idest,int isou)
    {
      vector_copy(levs[idest].conf,levs[isou].conf);
      levs[idest].nsmooth=levs[isou].nsmooth;
    }
    
    //evolve until nsmooth
    void update_level(int is,int nsmooth,int *dirs,int staple_min_dir)
    {
      while(levs[is].nsmooth<nsmooth)
	{
	  smooth_lx_conf_one_step(levs[is].conf,smoother,dirs,staple_min_dir);
	  levs[is].nsmooth++;
	}
    }
    
    //update the lowest conf using in chain the superior ones
    void update(int nsmooth,int *dirs,int staple_min_dir)
    {
      for(int is=1;is<nl();is++)
	{
	  int lev_ncur_smooth=levs[is].nsmooth;
	  int lev_ntarg_smooth=nsmooth%levs[is].epsb;
	  master_printf("Targeting at nsmooth=%d, current smooth at level %d: %d, this should be %d\n",nsmooth,is,lev_ncur_smooth,lev_ntarg_smooth);
	  
	  //if not on top, search the closest and extend it
	  if(lev_ncur_smooth!=lev_ntarg_smooth)
	    {
	      //find
	      int iclosest_smooth=find_closest_smaller_nsmooth(lev_ntarg_smooth);
	      master_printf("Closest level: %d, nsmooth: %d\n",iclosest_smooth,levs[iclosest_smooth].nsmooth);
	      
	      //copy and extend if needed
	      copy_level(is,iclosest_smooth);
	      update_level(is,nsmooth,dirs,staple_min_dir);
	    }
	}
    }
  };
}

#endif
