#ifndef GAUGE_SWEEP_H
#define GAUGE_SWEEP_H

#include "new_types/new_types_definitions.hpp"
#include "communicate/all_to_all.hpp"

namespace nissa
{
  struct gauge_sweep_t
  {
    bool path_inited,par_geom_inited;
    int comm_init_time;
    int nlinks_per_paths_site;
    int gpar;
    int max_cached_link;
    int max_sending_link;
    int *ilink_per_paths;
    int *nsite_per_box_dir_par;
    int *ivol_of_box_dir_par;
    all_to_all_comm_t *box_comm[16];
    su3 *buf_out,*buf_in;
    void(*add_paths_per_site_dir)(int *ilink_to_be_used,all_to_all_gathering_list_t &gat,int ivol,int mu);
    void init_box_dir_par_geometry(int ext_gpar,int(*par_comp)(coords ivol_coord,int dir));
    void init_paths(int ext_nlinks_per_paths_site,void(*ext_add_paths_per_site_dir)(int *ilink_to_be_used,all_to_all_gathering_list_t &gat,int ivol,int mu)
		    );
    void add_paths(all_to_all_gathering_list_t **gl);
    ~gauge_sweep_t();
    gauge_sweep_t();
  };
}

#endif
