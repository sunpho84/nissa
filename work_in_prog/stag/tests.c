#pragma once

//test the application of a backfield
int backfield_application_test()
{
  quad_su3 *eo_conf_comp[2]={nissa_malloc("conf_comp",loc_vol,quad_su3),eo_conf_comp[0]+loc_volh};
  quad_su3 *eo_conf_read_with_phases[2]={nissa_malloc("conf_read_with_phases",loc_vol,quad_su3),eo_conf_read_with_phases[0]+loc_volh};
  
  quad_u1 *u1b[2]={nissa_malloc("u1b",loc_vol,quad_u1),u1b[0]+loc_volh};
  
  ////////////////////////////////////////////
  
  //read the plain and phased conf
  read_ildg_conf_and_split_into_eo_parts(eo_conf_comp,"dat/conf_plain");
  read_ildg_conf_and_split_into_eo_parts(eo_conf_read_with_phases,"dat/conf_with_phases");
  
  //define phases
  init_backfield_to_id(u1b);
  add_stagphases_to_backfield(u1b);
  add_antiperiodic_bc_to_backfield(u1b);
  
  //add them
  add_backfield_to_conf(eo_conf_comp,u1b);
  
  //compute difference
  double glb_diff,loc_diff=0;
  for(int par=0;par<2;par++)
    for(int ivol=0;ivol<loc_volh;ivol++)
      for(int mu=0;mu<4;mu++)
	for(int ic1=0;ic1<3;ic1++)
	  for(int ic2=0;ic2<3;ic2++)
	    for(int ri=0;ri<2;ri++)
	      loc_diff+=pow(eo_conf_comp[par][ivol][mu][ic1][ic2][ri]-eo_conf_read_with_phases[par][ivol][mu][ic1][ic2][ri],2);
  MPI_Allreduce(&loc_diff,&glb_diff,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  
  glb_diff=sqrt(glb_diff/(2*3*3*4*loc_vol));
  
  master_printf("Difference between loaded and computed conf: %lg\n",glb_diff);
  
  ///////////////////////////////////////////
  
  nissa_free(u1b[0]);
  
  nissa_free(eo_conf_read_with_phases[0]);
  nissa_free(eo_conf_comp[0]);
  
  if(glb_diff<1.e-14) 
    {
      master_printf("Backfield application test passed!\n");
      return 1;
    }
  else
    {
      master_printf("Backfield application test not passed!\n");
      return 0;
    }
}

//test the application of D^2
int stD2ee_application_test()
{
  quad_su3 *eo_conf[2]={nissa_malloc("conf",loc_vol+loc_bord,quad_su3),eo_conf[0]+loc_volh+loc_bordh};
  quad_u1 *u1b[2]={nissa_malloc("u1b",loc_vol,quad_u1),u1b[0]+loc_volh};
  color *read_rnd=nissa_malloc("rnd",(loc_vol+loc_bord)/2,color);
  color *read_phi=nissa_malloc("read_phi",(loc_vol+loc_bord)/2,color);
  color *comp_phi=nissa_malloc("comp_phi",(loc_vol+loc_bord)/2,color);
    
  ////////////////////////////////////////////
  
  //read the plain and phased conf
  read_ildg_conf_and_split_into_eo_parts(eo_conf,"dat/conf_plain");
  
  //define phases
  init_backfield_to_id(u1b);
  add_stagphases_to_backfield(u1b);
  add_antiperiodic_bc_to_backfield(u1b);
  
  //add them
  add_backfield_to_conf(eo_conf,u1b);
  
  //read the rnd vector
  read_e_color(read_rnd,"dat/rnd_color_pfgen");
  
  //apply the dirac matr using read_phi as temp vector
  communicate_ev_color_borders(read_rnd);
  communicate_eo_gauge_borders(eo_conf[0],eo_conf[1]);
  apply_stD2ee(comp_phi,eo_conf,read_phi,0.025,read_rnd);
  
  //read the external phi
  read_e_color(read_phi,"dat/phi_color_md2ee_appl");
  
  //compute differnce between read and computed vector
  double glb_diff,loc_diff=0;
  for(int ivol=0;ivol<loc_volh;ivol++)
    for(int ic=0;ic<3;ic++)
      for(int ri=0;ri<2;ri++)
	loc_diff+=pow(comp_phi[ivol][ic][ri]-read_phi[ivol][ic][ri],2);
  MPI_Allreduce(&loc_diff,&glb_diff,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  glb_diff=sqrt(glb_diff/(loc_volh*3*2));
  
  nissa_free(comp_phi);
  nissa_free(read_phi);
  nissa_free(read_rnd);
  
  master_printf("Difference between loaded and computed vec: %lg\n",glb_diff);
  
  ///////////////////////////////////////////
  
  nissa_free(u1b[0]);
  
  nissa_free(eo_conf[0]);
  
  if(glb_diff<1.e-14) 
    {
      master_printf("stD2ee application test passed!\n");
      return 1;
    }
  else
    {
      master_printf("stD2ee application test not passed!\n");
      return 0;
    }
}

//test the application of (D2ee)^-(1/8)
int stD2ee_pow_minus_one_eighth_application_test()
{
  quad_su3 *eo_conf[2]={nissa_malloc("conf",loc_vol+loc_bord,quad_su3),eo_conf[0]+loc_volh+loc_bordh};
  quad_u1 *u1b[2]={nissa_malloc("u1b",loc_vol,quad_u1),u1b[0]+loc_volh};
  color *read_rnd=nissa_malloc("rnd",(loc_vol+loc_bord)/2,color);
  color *read_phi=nissa_malloc("read_phi",(loc_vol+loc_bord)/2,color);
  color *comp_phi=nissa_malloc("comp_phi",(loc_vol+loc_bord)/2,color);
    
  ////////////////////////////////////////////
  
  //read the plain and phased conf
  read_ildg_conf_and_split_into_eo_parts(eo_conf,"dat/conf_plain");
  
  //define phases
  init_backfield_to_id(u1b);
  add_stagphases_to_backfield(u1b);
  add_antiperiodic_bc_to_backfield(u1b);
  
  //read the rnd vector
  read_e_color(read_rnd,"dat/rnd_color_pfgen");
  
  //add them
  add_backfield_to_conf(eo_conf,u1b);
  
  //define the approximation of x^(-1/8)
  rat_approx appr;
  //load original approximation
  open_input("dat/rhmc4");
  read_double(&(appr.minimum));
  read_double(&(appr.maximum));
  int nterms;
  read_int(&nterms);
  rat_approx_create(&appr,nterms,"x^0.125");
  read_double(&appr.exp_power);
  read_double(&appr.cons);
  for(int iterm=0;iterm<nterms;iterm++)
    {
      read_double(&appr.poles[iterm]);
      read_double(&appr.weights[iterm]);
    }
  close_input();
  
  //print the approximation
  master_printf_rat_approx(&appr);
  
  //apply the whole rational approximation
  summ_src_and_all_inv_stD2ee_cgmm2s(comp_phi,read_rnd,eo_conf,&appr,1000000,1.e-24,1.e-24,0);
  
  //read the external phi
  read_e_color(read_phi,"dat/pow_minus_one_eighth_appl");
  
  //compute differnce between read and computed vector
  double glb_diff,loc_diff=0;
  for(int ivol=0;ivol<loc_volh;ivol++)
    for(int ic=0;ic<3;ic++)
      for(int ri=0;ri<2;ri++)
	loc_diff+=pow(comp_phi[ivol][ic][ri]-read_phi[ivol][ic][ri],2);
  MPI_Allreduce(&loc_diff,&glb_diff,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  glb_diff=sqrt(glb_diff/(loc_volh*3*2));
  
  nissa_free(comp_phi);
  nissa_free(read_phi);
  nissa_free(read_rnd);
  
  master_printf("Difference between loaded and computed vec: %lg\n",glb_diff);
  
  ///////////////////////////////////////////
  
  nissa_free(u1b[0]);
  
  nissa_free(eo_conf[0]);
  
  if(glb_diff<1.e-9) 
    {
      master_printf("stD2ee^(-1/8) application test passed!\n");
      return 1;
    }
  else
    {
      master_printf("stD2ee^(-1/8) application test not passed!\n");
      return 0;
    }
}

//test the determination of the maximal eigenvalue of D2ee
int stD2ee_max_eigenvalue_find_test()
{
  quad_su3 *eo_conf[2]={nissa_malloc("conf",loc_vol+loc_bord,quad_su3),eo_conf[0]+loc_volh+loc_bordh};
  quad_u1 *u1b[2]={nissa_malloc("u1b",loc_vol,quad_u1),u1b[0]+loc_volh};
    
  ////////////////////////////////////////////
  
  //read the plain and phased conf
  read_ildg_conf_and_split_into_eo_parts(eo_conf,"dat/conf_plain");
  
  //define phases
  init_backfield_to_id(u1b);
  add_stagphases_to_backfield(u1b);
  add_antiperiodic_bc_to_backfield(u1b);
  
  //add them
  add_backfield_to_conf(eo_conf,u1b);
  
  //find the maximal eigenvalue
  quark_content pars;
  pars.mass=0.025;
  double mc=max_eigenval(&pars,eo_conf,5000);
  
  ///////////////////////////////////////////
  
  nissa_free(u1b[0]);
  nissa_free(eo_conf[0]);
  
  /*
  if(glb_diff<1.e-9) 
    {
      master_printf("stD2ee^(-1/8) application test passed!\n");
      return 1;
    }
  else
    {
      master_printf("stD2ee^(-1/8) application test not passed!\n");
      return 0;
    }
  */
}



//generate pseudo-fermions using color vector generator
void generate_pseudo_fermions()
{
  /*
  color *temp=nissa_malloc("temp",loc_vol,color);  
  //generate the random field
  for(int ivol=0;ivol<loc_vol/2;ivol++)
    color_put_to_gauss(temp[ivol],&(loc_rnd_gen[ivol]),1);
  */
  
  /* 
  color *read_rnd=nissa_malloc("rnd",(loc_vol+loc_bord)/2,color);
  color *read_phi=nissa_malloc("read_phi",(loc_vol+loc_bord)/2,color);
  color *comp_phi=nissa_malloc("comp_phi",(loc_vol+loc_bord)/2,color);
  
  add_backfield_to_conf(eo_conf,u1b);
  communicate_eo_gauge_borders(eo_conf[0],eo_conf[1]);
  
  read_e_color(read_rnd,"rnd");
  communicate_ev_color_borders(read_rnd);
  
  summ_src_and_all_inv_stD2ee_cgmm2s(comp_phi,read_rnd,eo_conf,&(rat_exp_pfgen[0]),1000000,1.e-10,1.e-10,0);
  
  rem_backfield_to_conf(eo_conf,u1b);
  
  read_e_color(read_phi,"phie");
  
  double diff=0;
  for(int ivol=0;ivol<loc_vol/2;ivol++)
    for(int ic=0;ic<3;ic++)
      for(int ri=0;ri<2;ri++)
	{
	  master_printf("%lg %lg\n",comp_phi[ivol][ic][ri],read_phi[ivol][ic][ri]);
	  diff+=pow(comp_phi[ivol][ic][ri]-read_phi[ivol][ic][ri],2);
	}
  diff/=loc_vol/2*2*3;

  master_printf("%lg\n",sqrt(diff));
  
  nissa_free(comp_phi);
  nissa_free(read_phi);
  nissa_free(read_rnd);
 */
  //  nissa_free(temp);
}
