
//multiply an su2 matrix and an su3 and assign to last
void unsafe_su2_prodassign_su3(int ic1,int ic2,double a0,double a1,double a2,double a3,su3 in)
{
  complex c11={ a0, a3};
  complex c22={ a0,-a3};
  complex c12={ a2, a1};
  complex c21={-a2, a1};
  
  complex row1[3];
  complex row2[3];
  
  //create the two new rows of the matrix
  for(int ic=0;ic<3;ic++)
    {
      //first row
      safe_complex_prod    (row1[ic],c11,in[ic1][ic]);
      complex_summ_the_prod(row1[ic],c12,in[ic2][ic]);
      //second row
      safe_complex_prod    (row2[ic],c21,in[ic1][ic]);
      complex_summ_the_prod(row2[ic],c22,in[ic2][ic]);
    }
    
  //change the two lines in the matrix
  for(int ic=0;ic<3;ic++)
    {
      complex_copy(in[ic1][ic],row1[ic]);
      complex_copy(in[ic2][ic],row2[ic]);
    }
}

//return a single link after the heatbath procedure
void find_heatbath_link(su3 u,quad_su3 **eo_conf,int par,int ieo,int mu,double beta,int hb_steps)
{
  int ilx=loclx_of_loceo[par][ieo];
  rnd_gen *gen=&(loc_rnd_gen[ilx]);
  
  //subgroups indices
  int ic[3][2]={{0,1},{0,2},{1,2}};
  
  //compute the staple
  su3 staple;
  compute_point_staples_eo_conf_single_dir(staple,eo_conf,ilx,mu);
  
  //compute the original contribution to the action due to the given link 
  su3 uprod;
  unsafe_su3_prod_su3_dag(uprod,eo_conf[par][ieo][mu],staple);
  
  for(int hb=0;hb<hb_steps;hb++)
    {
      //extract the subgroup to use                                                                                                                                             
      int igr=(int)rnd_get_unif(gen,0,3);
      
      //color indeces defining subgroup
      int ic1=ic[igr][0];
      int ic2=ic[igr][1];
      
      complex c1,c2;
      complex_summ_conj2(c1,uprod[ic1][ic1],uprod[ic2][ic2]);
      complex_subt_conj2(c2,uprod[ic1][ic2],uprod[ic2][ic1]);
      
      double smod=1/sqrt(c1[0]*c1[0]+c1[1]*c1[1]+c2[0]*c2[0]+c2[1]*c2[1]);                                                                                                      
      double omega_f=beta/(3*smod);
      
      double r0=c1[RE]*smod;
      double r1=c2[IM]*smod;
      double r2=c2[RE]*smod;
      double r3=c1[IM]*smod;
      
      double z_norm=exp(-2*omega_f);
      omega_f=1/omega_f;
      
      double temp_f,z_f,a0;
      do
	{
	  double z_temp = (z_norm-1)*rnd_get_unif(gen,0,1)+1;
	  a0     = 1+omega_f*log(z_temp);
	  z_f    = 1-a0*a0;
	  temp_f = sqr(rnd_get_unif(gen,0,1))-z_f;
	}
      while(temp_f>0);
      
      double x_rat=sqrt(z_f);
      
      //generate an su2 matrix
      double fi=rnd_get_unif(gen,0,2*M_PI);
      double cteta=rnd_get_unif(gen,-1,1);
      double steta=sqrt(1-cteta*cteta);
      
      double a1=steta*cos(fi)*x_rat;
      double a2=steta*sin(fi)*x_rat;
      double a3=cteta*x_rat;
      
      double x0 = a0*r0 + a1*r1 + a2*r2 + a3*r3;
      double x1 = r0*a1 - a0*r1 + a2*r3 - r2*a3;
      double x2 = r0*a2 - a0*r2 + a3*r1 - r3*a1;
      double x3 = r0*a3 - a0*r3 + a1*r2 - r1*a2;

      unsafe_su2_prodassign_su3(ic1,ic2,x0,x1,x2,x3,uprod);
      su3_copy(u,eo_conf[par][ieo][mu]);
      unsafe_su2_prodassign_su3(ic1,ic2,x0,x1,x2,x3,u);
    }
}

//heatbath algorithm for the quenched simulation case
void heatbath(quad_su3 **eo_conf,double beta,int hb_steps)
{
  //loop first on parity and then on directions                                                                                                                           
  for(int mu=0;mu<4;mu++)
    for(int par=0;par<2;par++)
      {
        nissa_loc_volh_loop(ieo)
	{
	  //find the transformation                                                                                                                                           
	  su3 u;
	  find_heatbath_link(u,eo_conf,par,ieo,mu,beta,hb_steps);
	  
	  //apply it
	  su3_copy(eo_conf[par][ieo][mu],u);
	}
	
        //now set the borders invalid: since we split conf in e/o, only now needed                                                                                              
        set_borders_invalid(eo_conf[par]);
      }
}
