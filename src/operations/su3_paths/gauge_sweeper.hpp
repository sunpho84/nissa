#ifndef _GAUGE_SWEEPER_HPP
#define _GAUGE_SWEEPER_HPP

#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "base/vectors.hpp"
#include "communicate/all_to_all.hpp"
#include "hmc/gauge/gluonic_action.hpp"
#include "new_types/su3.hpp"
#include "new_types/su3_op.hpp"
#include "threads/threads.hpp"

#ifndef EXTERN
 #define EXTERN extern
#endif

namespace nissa
{
  //sweep a configuration, possibly using subboxes, each divided in checkboard so to avoid communication problem
  struct gauge_sweeper_t
  {
    //flags
    bool staples_inited,par_geom_inited,packing_inited;
    
    //coefficient for rectangles
    //this is fixed when asking a sweeper
    double C1;
    
    //benchmarks and checks
    double comm_init_time,comp_time,comm_time;
    int max_cached_link,max_sending_link;
    
    //store action parameters
    int nlinks_per_staples_of_link,gpar;
    
    //alternative ways to compute
    int *ilink_per_staples;
    int *packing_link_source_dest;//std::map<int,std::vector<int> > *packing_index;
    su3 *packing_link_buf;
    //geometry
    int *nsite_per_box_dir_par;
    int *ivol_of_box_dir_par;
    
    //communicators
    all_to_all_comm_t *box_comm[16];
    su3 *buf_out,*buf_in;
    
    ///////////////////////////////// methods ///////////////////////
    
    //routine used to add paths (pointer to external function is stored here for thread commodity used)
    void(*add_staples_per_link)(int *ilink_to_be_used,all_to_all_gathering_list_t &gat,int ivol,int mu);
    void init_staples(int ext_nlinks_per_staples_of_link,void(*ext_add_staples_per_link)
                      (int *ilink_to_be_used,all_to_all_gathering_list_t &gat,int ivol,int mu),
                      void (*ext_compute_staples)(su3 staples,su3 *links,int *ilinks,double C1));
    void add_staples_required_links(all_to_all_gathering_list_t **gl);
    
    //find the order in which to scan the links to compute the staple sequentially
    void find_packing_index(void (*ext_compute_staples_packed)(su3 staples,su3 *links,double C1));
    void pack_links(quad_su3 *conf,int ibase,int nbox_dir_par);
    
    //routine computing staples
    void (*compute_staples)(su3 staples,su3 *links,int *ilinks,double C1);
    void (*compute_staples_packed)(su3 staples,su3 *links,double C1);
    
    //inits the parity checkboard according to an external parity
    void init_box_dir_par_geometry(int ext_gpar,int(*par_comp)(const Coords& ivol_coord,const int& dir));
    
    //sweep the conf
    void sweep_conf(LxField<quad_su3>& conf,void (*update_fun)(su3 out,su3 staples,int ivol,int mu,void *pars),void *pars)
    {
      CRASH("reimplement");
      // MANDATORY_PARALLEL;
      
      // int ibase=0;
      // for(int ibox=0;ibox<(1<<NDIM);ibox++)
      // 	{
      // 	  //communicate needed links
      // 	  if(IS_MASTER_THREAD) comm_time-=take_time();
      // 	  box_comm[ibox]->communicate(conf,conf,sizeof(su3),NULL,NULL,ibox+100);
      // 	  if(IS_MASTER_THREAD)
      // 	    {
      // 	      comm_time+=take_time();
      // 	      comp_time-=take_time();
      // 	    }
      // 	  for(int dir=0;dir<NDIM;dir++)
      // 	    for(int par=0;par<gpar;par++)
      // 	      {
      // 		int box_dir_par_size=nsite_per_box_dir_par[par+gpar*(dir+NDIM*ibox)];
		
      // 		//pack
      // 		if(packing_inited) pack_links(conf,ibase,box_dir_par_size);
		
      // 		//scan the whole box
      // 		NISSA_PARALLEL_LOOP(ibox_dir_par,ibase,ibase+box_dir_par_size)
      // 		  {
      // 		    //compute the staples
      // 		    su3 staples;
		    
      // 		    if(packing_inited)
      // 		      compute_staples_packed(staples,packing_link_buf+(ibox_dir_par-ibase)*nlinks_per_staples_of_link,C1);
      // 		    else
      // 		      compute_staples(staples,(su3*)conf,ilink_per_staples+nlinks_per_staples_of_link*ibox_dir_par,C1);
		    
      // 		    //find new link
      // 		    int ivol=ivol_of_box_dir_par[ibox_dir_par];
      // 		    update_fun(conf[ivol][dir],staples,ivol,dir,pars);
      // 		  }
      // 		NISSA_PARALLEL_LOOP_END;
      // 		THREAD_BARRIER();
		
      // 		//increment the box-dir-par subset
      // 		ibase+=box_dir_par_size;
      // 	      }
      // 	  if(IS_MASTER_THREAD) comp_time+=take_time();
      // 	}
      
      // set_borders_invalid(conf);
    }
  
    //checkers
    void check_hit_in_the_exact_order();
    void check_hit_exactly_once();
    
    //constructor, destructors
    ~gauge_sweeper_t();
    gauge_sweeper_t();
  };
  
  void compute_Symanzik_staples_packed(su3 staples,su3 *links,double C1);
  void init_Symanzik_sweeper();
  void init_Wilson_sweeper();
  void init_sweeper(gauge_action_name_t);
  
  EXTERN gauge_sweeper_t *Symanzik_sweeper,*Wilson_sweeper;
  
  inline gauge_sweeper_t *get_sweeper(gauge_action_name_t gauge_action_name)
  {
    gauge_sweeper_t * sweeper=NULL;
    
    switch(gauge_action_name)
      {
      case WILSON_GAUGE_ACTION:
	init_Wilson_sweeper();
	Wilson_sweeper->C1=0;
	sweeper=Wilson_sweeper;
	break;
      case TLSYM_GAUGE_ACTION:
	init_Symanzik_sweeper();
	Symanzik_sweeper->C1=C1_TLSYM;
	sweeper=Symanzik_sweeper;
      case IWASAKI_GAUGE_ACTION:
	init_Symanzik_sweeper();
	Symanzik_sweeper->C1=C1_IWASAKI;
	sweeper=Symanzik_sweeper;
	break;
      case UNSPEC_GAUGE_ACTION:
	CRASH("unspecified action");
      default:
	CRASH("not implemented action");
	sweeper=NULL;
	break;
      }
    
    return sweeper;
  }
}

#undef EXTERN

#endif
