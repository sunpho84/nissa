#include "nissa.h"

int main(int narg,char **arg)
{
  //basic mpi initialization
  init_nissa();
  
  //open input file
  open_input("input");

  //set sizes
  int L,T;
  read_str_int("L",&L);
  read_str_int("T",&T);
  
  //Init the MPI grid 
  init_grid(T,L);
  
  close_input();
  
  ///////////////////////////////////////
  
  quad_su3 *lx_conf=nissa_malloc("lx_conf",loc_vol+loc_bord+loc_edge,quad_su3);
  quad_su3 *eo_conf[2]={nissa_malloc("ev_conf",loc_volh+loc_bordh+loc_edgeh,quad_su3),nissa_malloc("od_conf",loc_volh+loc_bordh+loc_edgeh,quad_su3)};

  read_ildg_gauge_conf(lx_conf,"dat/conf_plain");
  split_lx_conf_into_eo_parts(eo_conf,lx_conf);
  
  ///////////////////////////////////////
  
  communicate_lx_quad_su3_borders(lx_conf);
  communicate_eo_quad_su3_borders(eo_conf);
  
  for(int ibulk=0;ibulk<loc_vol;ibulk++)
    for(int vers=0;vers<2;vers++)
    for(int mu=0;mu<4;mu++)
      {
	int ibord=loclx_neigh[vers][ibulk][mu];
	if(ibord>=loc_vol)
	  {
	    int par=loclx_parity[ibord];
	    int ieo=loceo_of_loclx[ibord];
	    
	    if(memcmp(lx_conf[ibord],eo_conf[par][ieo],sizeof(quad_su3))) crash("error %d",ibord);
	  }
      }
  
  ///////////////////////////////////////
  
  communicate_lx_quad_su3_edges(lx_conf);
  communicate_eo_quad_su3_edges(eo_conf);
  
  for(int ibord=0;ibord<loc_vol+loc_bord;ibord++)
    for(int vers=0;vers<2;vers++)
    for(int mu=0;mu<4;mu++)
      {
	int iedge=loclx_neigh[vers][ibord][mu];
	if(iedge>=loc_vol+loc_bord)
	  {
	    int par=loclx_parity[iedge];
	    int ieo=loceo_of_loclx[iedge];
	    
	    if(memcmp(lx_conf[iedge],eo_conf[par][ieo],sizeof(quad_su3))) crash("error %d",iedge);
	  }
      }
  
  ///////////////////////////////////////
  
  nissa_free(lx_conf);
  nissa_free(eo_conf[0]);
  nissa_free(eo_conf[1]);
  
  close_nissa();
  
  return 0;
}
