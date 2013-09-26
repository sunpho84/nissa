#pragma once

int test_Q2tm_inversion()
{
  //load the well known source
  master_printf("\nLoading conf\n");
  quad_su3 *conf=nissa_malloc("conf",loc_vol+bord_vol+edge_vol,quad_su3);
  read_ildg_gauge_conf(conf,"../../data/L4T8conf");
  
  //generate the classic source
  master_printf("Generating source\n");
  spincolor *source=nissa_malloc("source",loc_vol+bord_vol,spincolor);
  if(nissa_loc_rnd_gen_inited) stop_loc_rnd_gen();
  start_loc_rnd_gen(2342);
  generate_undiluted_source(source,RND_Z4,-1);
  
  //invert Q2
  double kappa=0.177000;
  double mu=0.50;
  double prec=1.e-25;
  spincolor *inver=nissa_malloc("inver",loc_vol+bord_vol,spincolor);
  inv_Q2_cg(inver,source,NULL,conf,kappa,mu,1000000,5,prec);
  
  //now compare with saved data
  master_printf("Reading saved spincolor\n");
  spincolor *comp=nissa_malloc("comp",loc_vol,spincolor);
  read_spincolor(comp,"../../data/Q2tm_inv");
  
  //compare the weighted norm
  double loc_weighted_norm=0,weighted_norm;
  NISSA_LOC_VOL_LOOP(ivol)
    for(int id=0;id<4;id++)
      for(int ic=0;ic<3;ic++)
	for(int ri=0;ri<2;ri++)
	  loc_weighted_norm+=sqr((comp[ivol][id][ic][ri]-inver[ivol][id][ic][ri])/
				 (comp[ivol][id][ic][ri]+inver[ivol][id][ic][ri]));
  MPI_Allreduce(&loc_weighted_norm,&weighted_norm,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  weighted_norm=sqrt(weighted_norm/glb_vol/12);
  
  double tolerance=1.e-13;
  master_printf("Difference: %lg, tolerance: %lg\n",weighted_norm,tolerance);
  
  nissa_free(comp);
  nissa_free(inver);
  nissa_free(source);
  nissa_free(conf);
  
  return weighted_norm<=tolerance;
}
