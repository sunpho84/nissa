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
  //holds the data to recursive Wflow, generalizing from 1302.5246
  struct recursive_Wflower_t
  {
    struct Wflow_lev_t
    {
      quad_su3 *conf; //stored configuration
      int nWflow; //current Wflow in terms of the fundamental
      int units; //length of the Wflowing in units of the fundamental, for this level
      int off() //offset of the level
      {return units/2;}
    };
    
    const Wflow_pars_t &Wflower;
    
    //level of Wflowed conf
    std::vector<Wflow_lev_t> levs;
    
    //number of levels
    int nl()
    {return levs.size();}
    
    //initialize with given flower and conf
    recursive_Wflower_t(const Wflow_pars_t &Wflower,quad_su3 *ori_conf) : Wflower(Wflower)
    {
      int ns=Wflower.nrecu;
      master_printf("ns: %d\n",ns);
      
      //take the number of Wflow levels
      int m=Wflower.nflows;
      master_printf("m: %d\n",m);
      //check
      if(ns<0 or ns>m) crash("ns=%d must be in the range [%d,%d]",0,m);
      
      //set the optimal number of blocks
      double nb=pow(m,1.0/ns);
      //verbosity_lv3_
	master_printf("nb: %lg\n",nb);
      
      //resize
      levs.resize(ns+1);
      for(int i=0;i<ns;i++) levs[i].conf=nissa_malloc("conf",(locVol+bord_vol+edge_vol).nastyConvert(),quad_su3);
      
      //set units of Wflow
      for(int il=0;il<nl();il++)
	{
	  levs[il].units=std::max(1.0,m/pow(nb,il)+1e-10);
	  // verbosity_lv3_
	    master_printf("unit %d: %d\n",il,levs[il].units);
	}
      if(levs.back().units!=1) crash("units of the lowest level has to be 1, instead it is %d",levs.back().units);
      
      bind_conf(ori_conf);
    }
    
    //destroyer
    ~recursive_Wflower_t()
    {for(int is=0;is<nl()-1;is++) nissa_free(levs[is].conf);}
    
    //set the external conf
    void bind_conf(quad_su3 *ori_conf)
    {
      master_printf("binding\n");
      //set the pointer of the lowest level to the external (to be evolved) conf
      levs.back().conf=ori_conf;
      //set units of toppest to a very large number so it's never evolved
      levs.front().units=200000;
      for(int i=0;i<nl();i++)
	{
	  master_printf("copying i %d\n",i);
	  vector_copy(levs[i].conf,ori_conf);
	  levs[i].nWflow=0;
	}
    }
    
    //find the level which has the closest smaller nWflow
    int find_closest_smaller_nWflow(int nWflow)
    {
      //verbosity_lv3_
	master_printf(" searching the level closest to %d\n",nWflow);
      
      int iclosest_Wflow=0;
      for(int is=0;is<nl();is++)
	{
	  int nclosest_Wflow=levs[iclosest_Wflow].nWflow;
	  int lev_ncur_Wflow=levs[is].nWflow;
	  //verbosity_lv3_
	    master_printf(" level %d: %d, current closest: %d, nclosest_Wflow: %d\n",is,lev_ncur_Wflow,iclosest_Wflow,nclosest_Wflow);
	  if(lev_ncur_Wflow<=nWflow and lev_ncur_Wflow>nclosest_Wflow) iclosest_Wflow=is;
	}
      
      return iclosest_Wflow;
    }
    
    //copy level
    void copy_level(int idest,int isou)
    {
      vector_copy(levs[idest].conf,levs[isou].conf);
      levs[idest].nWflow=levs[isou].nWflow;
    }
    
    //evolve until nWflow
    int update_level(int is,int nWflow,bool *dirs=all_dirs)
    {
      int nevol=0;
      
      //verbosity_lv3_
	master_printf("levs[%d].nWflow: %d nWflow: %d\n",is,levs[is].nWflow,nWflow);
      while(levs[is].nWflow<nWflow)
	{
	  Wflow_lx_conf(levs[is].conf,Wflower.dt,dirs);
	  levs[is].nWflow++;
	  nevol++;
	}
      
      return nevol;
    }
    
    //return the smallest viable nWflow rounded to the level
    int get_smallest_nWflow_rounded(int nWflow,int is)
    {
      int units=levs[is].units;
      int off=levs[is].off();
      int targ=((nWflow+units-off)/units-1)*units+off;
      if(targ<0) targ=0;
      return targ;
    }
    
    //update the lowest conf using in chain the superior ones
    int update(int nWflow,bool *dirs=all_dirs,int staple_min_dir=0)
    {
      int nevol=0;
      
      //store the path
      std::vector<std::pair<int,int> > order_path;
      
      for(int is=1;is<nl();is++)
	{
	  int lev_ncur_Wflow=levs[is].nWflow;
	  int lev_ntarg_Wflow=get_smallest_nWflow_rounded(nWflow,is);
	  int units=levs[is].units;
	  //verbosity_lv3_
	    master_printf("Targeting at nWflow=%d, current Wflow at level %d units of %d: %d, this should be %d\n",nWflow,is,units,lev_ncur_Wflow,lev_ntarg_Wflow);
	    
	    //store in the path the current level, using the target nWflow for current level as a key
	    if(lev_ncur_Wflow!=lev_ntarg_Wflow and lev_ntarg_Wflow>=0 and lev_ntarg_Wflow<=nWflow)
	      order_path.push_back(std::make_pair(lev_ntarg_Wflow,is));
	}
      
      //sort the path
      std::sort(order_path.begin(),order_path.end());
      for(int ip=0;ip<(int)order_path.size();ip++) verbosity_lv3_master_printf("ip=%d targ=%d is=%d\n",ip,order_path[ip].first,order_path[ip].second);
      
      //if not on correct number of Wflow, search for the closest and extend it
      for(int ip=0;ip<(int)order_path.size();ip++)
	{
	  int lev_ntarg_Wflow=order_path[ip].first;
	  int il=order_path[ip].second;
	  
	  int iclosest_Wflow=find_closest_smaller_nWflow(lev_ntarg_Wflow);
	  //verbosity_lv3_
	  master_printf("Targetting %d for lev %d, closest level: %d, nWflow: %d\n",lev_ntarg_Wflow,il,iclosest_Wflow,levs[iclosest_Wflow].nWflow);
	  
	  if(iclosest_Wflow!=il)
	    {
	      master_printf("Copying from %d\n",iclosest_Wflow);
       	      copy_level(il,iclosest_Wflow);
	    }
	  
	  //extend if needed
	  if(levs[il].nWflow!=lev_ntarg_Wflow) nevol+=update_level(il,lev_ntarg_Wflow,dirs);
	  else master_printf("Level %d alright\n",il);
       	}
      
      return nevol;
    }
  };
}

#endif
