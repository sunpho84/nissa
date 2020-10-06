#pragma once

//initialize e/o geometry
int test_eo_geometry()
{
  const double tolerance=1.e-14;
  
  //set_eo_geometry();
  
  //allocate lx and eo conf
  quad_su3 *lx_conf=nissa_malloc("lx_conf",loc_vol+bord_vol+edge_vol,quad_su3);
  quad_su3 *eo_conf[2];
  eo_conf[0]=nissa_malloc("eo_conf",loc_vol+bord_vol+edge_vol,quad_su3);
  eo_conf[1]=eo_conf[0]+(loc_vol+bord_vol+edge_vol)/2;
  quad_su3 *lx_conf_reco=nissa_malloc("lx_conf_reco",loc_vol+bord_vol+edge_vol,quad_su3);
  
  //read conf, compute plaquette, print it
  read_ildg_gauge_conf(lx_conf,"../../data/L4T8conf");
  communicate_lx_quad_su3_borders(lx_conf);
  communicate_lx_gauge_edges(lx_conf);
  
  split_lx_vector_into_eo_parts((char*)(eo_conf[0]),(char*)(eo_conf[1]),(char*)lx_conf,loc_vol,sizeof(quad_su3));
  paste_eo_parts_into_lx_vector((char*)(lx_conf_reco),(char*)(eo_conf[0]),(char*)eo_conf[1],loc_vol+bord_vol+edge_vol,sizeof(quad_su3));
  
  double plaq_lx=global_plaquette_lx_conf(lx_conf);
  double plaq_eo=global_plaquette_eo_conf(eo_conf[0],eo_conf[1]);
  communicate_ev_and_od_quad_su3_borders(eo_conf[0],eo_conf[1]);
  double plaq_eo_eo_comm=global_plaquette_eo_conf(eo_conf[0],eo_conf[1]);
  double plaq_reco=global_plaquette_lx_conf(lx_conf_reco);
  master_printf("lx plaquette: %16.16lg, eo comp %16.16lg, eo comm comp %16.16lg, reconstructed: %16.16lg\n",plaq_lx,plaq_eo,plaq_eo_eo_comm,plaq_reco);
  
  int passed1=fabs(plaq_lx-plaq_eo)<=tolerance;
  if(!passed1) master_printf("lx and eo plaquette do not agree!\n");
  int passed2=fabs(plaq_lx-plaq_reco)<=tolerance;
  if(!passed2) master_printf("lx direct and eo split-reconstruct do not agree!\n");
  int passed3=fabs(plaq_eo_eo_comm-plaq_lx)<=tolerance;
  if(!passed3) master_printf("lx direct and eo communicate border do not agree!\n");

  nissa_free(lx_conf_reco);
  nissa_free(eo_conf[0]);
  nissa_free(lx_conf);
  
  unset_eo_geometry();
  
  return passed1 && passed2 && passed3;
}
