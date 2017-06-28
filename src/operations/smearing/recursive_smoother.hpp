#ifndef _RECURSIVE_SMOOTHER_HPP
#define _RECURSIVE_SMOOTHER_HPP

#include "base/vectors.hpp"
#include "geometry/geometry_lx.hpp"
#include "new_types/su3.hpp"
#include "operations/smearing/smooth.hpp"

#include <utility>
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
      int units; //length of the smoothing in units of the fundamental, for this level
      int off() //offset of the level
      {return units/2;}
    };
    
    smooth_pars_t &smoother; //smoother embedded
    
    //level of smoothed conf
    std::vector<smooth_lev_t> levs;
    
    //number of levels
    int nl()
    {return levs.size();}
    
    recursive_smoother_t(int ns,smooth_pars_t &smoother) : smoother(smoother)
    {
      //take the number of smooth levels
      int m=smoother.nsmooth();
      
      //check
      if(ns<0 or ns>m) crash("ns=%d must be in the range [%d,%d]",0,m);
      
      //set the optimal number of blocks
      double nb=pow(m,1.0/(ns+1.0));
      verbosity_lv3_master_printf("nb: %lg\n",nb);
      
      //resize
      levs.resize(ns+1);
      for(int i=0;i<nl();i++) levs[i].conf=nissa_malloc("conf",loc_vol+bord_vol,quad_su3);
      
      //set units of smooth
      for(int il=0;il<nl();il++)
	{
	  levs[il].units=std::max(1.0,m/pow(nb,il+1)+1e-10);
	  verbosity_lv3_master_printf("%d %.16lg,%d\n",il,std::max(1.0,m/pow(nb,il+1)),levs[il].units);
	}
    }
    
    ~recursive_smoother_t()
    {for(int is=0;is<nl();is++) nissa_free(levs[is].conf);}
    
    //set the external conf
    void set_conf(quad_su3 *ori_conf)
    {
      for(int i=0;i<nl();i++)
	{
	  vector_copy(levs[i].conf,ori_conf);
	  levs[i].nsmooth=-1;
	}
    }
    
    //find the level which has the closest smaller nsmooth
    int find_closest_smaller_nsmooth(int nsmooth)
    {
      //verbosity_lv3_
	master_printf(" searching the level closest to %d\n",nsmooth);
      
      int iclosest_smooth=0;
      for(int is=1;is<nl();is++)
	{
	  int nclosest_smooth=levs[iclosest_smooth].nsmooth;
	  int lev_ncur_smooth=levs[is].nsmooth;
	  //verbosity_lv3_
	    master_printf(" level %d: %d, current closest: %d, nclosest_smooth: %d\n",is,lev_ncur_smooth,iclosest_smooth,nclosest_smooth);
	  if(lev_ncur_smooth<=nsmooth and lev_ncur_smooth>nclosest_smooth) iclosest_smooth=is;
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
    void update_level(int is,int nsmooth,int *dirs=all_dirs,int staple_min_dir=0)
    {
      //verbosity_lv3_
	master_printf("levs[%d].nsmooth: %d nsmooth: %d\n",is,levs[is].nsmooth,nsmooth);
      while(levs[is].nsmooth<nsmooth)
	{
	  smooth_lx_conf_one_step(levs[is].conf,smoother,dirs,staple_min_dir);
	  levs[is].nsmooth++;
	}
    }
    
    //return the smallest viable nsmooth rounded to the level
    int get_smallest_nsmooth_rounded(int nsmooth,int is)
    {
      int units=levs[is].units;
      int off=levs[is].off();
      return ((nsmooth+units-off)/units-1)*units+off;
    }
    
    //update the lowest conf using in chain the superior ones
    quad_su3* update(int nsmooth,int *dirs=all_dirs,int staple_min_dir=0)
    {
      //store the path
      std::vector<std::pair<int,int> > order_path;
      
      for(int is=1;is<nl();is++)
	{
	  int lev_ncur_smooth=levs[is].nsmooth;
	  int lev_ntarg_smooth=get_smallest_nsmooth_rounded(nsmooth,is);
	  int units=levs[is].units;
	  //verbosity_lv3_
	    master_printf("Targeting at nsmooth=%d, current smooth at level %d units of %d: %d, this should be %d\n",nsmooth,is,units,lev_ncur_smooth,lev_ntarg_smooth);
	    
	    //store in the path the current level, using the target nsmooth for current level as a key
	    if(lev_ntarg_smooth>0 and lev_ntarg_smooth<=nsmooth)
	      order_path.push_back(std::make_pair(lev_ntarg_smooth,is));
	}
      
      //sort the path
      std::sort(order_path.begin(),order_path.end());
      for(int ip=0;ip<(int)order_path.size();ip++) verbosity_lv3_master_printf("ip=%d targ=%d is=%d\n",ip,order_path[ip].first,order_path[ip].second);
      
      // //if not on correct number of smooth, search for the closest and extend it
      // if(lev_ncur_smooth!=lev_ntarg_smooth and lev_ntarg_smooth<nsmooth)
      // 	    {
      // 	      //find
      // 	      int iclosest_smooth=find_closest_smaller_nsmooth(lev_ntarg_smooth);
      // 	      //verbosity_lv3_
      // 		master_printf("Closest level: %d, nsmooth: %d\n",iclosest_smooth,levs[iclosest_smooth].nsmooth);
	      
      // 	      //copy and extend if needed
      // 	      copy_level(is,iclosest_smooth);
      // 	      update_level(is,nsmooth,dirs,staple_min_dir);
      // 	    }
      // 	}
      
      return levs.back().conf;
    }
  };
}

#endif
