#include <mpi.h>
#include <lemon.h>

#include "appretto.h"

const int debug3=0;

double kappa;
double mu;
double muh;

void assign_add_mul_r_ee(spincolor *r,spincolor *s,double c)
{
  for(int X=0;X<loc_vol;X++)
    if(loclx_parity[X]==0)
      for(int id=0;id<4;id++)
	for(int ic=0;ic<3;ic++)
	  for(int ri=0;ri<2;ri++)
	    r[X][id][ic][ri]+=s[X][id][ic][ri]*c;
}

void assign_mul_one_pm_imu_inv_oo_or_ee(spincolor *psi,spincolor *source,double _sign,int parity)
{
  double sign=-1.; 
  if(_sign < 0.){
    sign = 1.; 
  }
  
  double nrm = 1/(1+muh*muh);
  complex z={nrm,sign*nrm*muh};
  complex w={z[0],-z[1]};
  spincolor temp;

  for(int X=0;X<loc_vol;X++)
    if(loclx_parity[X]==parity)
      {
	for(int ic=0;ic<3;ic++)
	  {
	    unsafe_complex_prod(temp[0][ic],source[X][0][ic],z);
	    unsafe_complex_prod(temp[1][ic],source[X][1][ic],z);
	    unsafe_complex_prod(temp[2][ic],source[X][2][ic],w);
	    unsafe_complex_prod(temp[3][ic],source[X][3][ic],w);
	  }
	memcpy(psi[X],temp,sizeof(spincolor));
      }
}

void assign_mul_one_pm_imu_inv_oo(spincolor *psi,spincolor *source,double sign)
{assign_mul_one_pm_imu_inv_oo_or_ee(psi,source,sign,1);}
void assign_mul_one_pm_imu_inv_ee(spincolor *psi,spincolor *source,double sign)
{assign_mul_one_pm_imu_inv_oo_or_ee(psi,source,sign,0);}

void assign_mul_one_pm_imu_sub_mul_gamma5_oo(spincolor *t,spincolor *r,spincolor *s,double _sign)
{
  double sign=1.; 
  if(_sign < 0.){
    sign = -1.; 
  }
  
  complex z={1,sign*muh};
  complex w={1,-sign*muh};
  spin phi;

  for(int X=0;X<loc_vol;X++)
    if(loclx_parity[X]==1)
      for(int ic=0;ic<3;ic++)
	{
	  unsafe_complex_prod(phi[0],z,r[X][0][ic]);
	  unsafe_complex_prod(phi[1],z,r[X][1][ic]);
	  unsafe_complex_prod(phi[2],w,r[X][2][ic]);
	  unsafe_complex_prod(phi[3],w,r[X][3][ic]);
	  
	  complex_subt(t[X][0][ic],phi[0],s[X][0][ic]);
	  complex_subt(t[X][1][ic],phi[1],s[X][1][ic]);
	  complex_subt(t[X][2][ic],s[X][2][ic],phi[2]);
	  complex_subt(t[X][3][ic],s[X][3][ic],phi[3]);
	}
}

void assign_mul_one_pm_imu_oo_or_ee(spincolor *psi,spincolor *source,double _sign,int parity)
{
  double sign=1.; 
  if(_sign < 0.){
    sign = -1.; 
  }
  
  complex z={1,sign*muh};
  complex w={z[0],-z[1]};
  spincolor temp;
  
  for(int X=0;X<loc_vol;X++)
    if(loclx_parity[X]==parity)
      {
	for(int ic=0;ic<3;ic++)
	  {
	    unsafe_complex_prod(temp[0][ic],source[X][0][ic],z);
	    unsafe_complex_prod(temp[1][ic],source[X][1][ic],z);
	    unsafe_complex_prod(temp[2][ic],source[X][2][ic],w);
	    unsafe_complex_prod(temp[3][ic],source[X][3][ic],w);
	  }
	memcpy(psi[X],temp,sizeof(spincolor));
      }
}

void assign_mul_one_pm_imu_oo(spincolor *psi,spincolor *source,double sign)
{assign_mul_one_pm_imu_inv_oo_or_ee(psi,source,sign,1);}
void assign_mul_one_pm_imu_ee(spincolor *psi,spincolor *source,double sign)
{assign_mul_one_pm_imu_inv_oo_or_ee(psi,source,sign,0);}

void assign_mul_add_r_oo_or_ee(spincolor *r,double c,spincolor *s,int parity)
{
  for(int X=0;X<loc_vol;X++)
    if(loclx_parity[X]==parity)
      for(int id=0;id<4;id++)
	for(int ic=0;ic<3;ic++)
	  for(int ri=0;ri<2;ri++)
	    r[X][id][ic][ri]=c*r[X][id][ic][ri]+s[X][id][ic][ri];
}

void assign_mul_add_r_oo(spincolor *r,double c,spincolor *s)
{assign_mul_add_r_oo_or_ee(r,c,s,1);}
void assign_mul_add_r_ee(spincolor *r,double c,spincolor *s)
{assign_mul_add_r_oo_or_ee(r,c,s,0);}

//apply D_eo or oe
void apply_Doe_or_Deo(spincolor *out,spincolor *in,quad_su3 *conf,int parity)
{
  int Xup,Xdw;

  if(debug3) printf(" m2Doe_or_m2Deo parity: %d\n",parity);
  for(int X=0;X<loc_vol;X++)
    if(loclx_parity[X]!=parity) //loop is on sink
      {
	color temp_c0,temp_c1,temp_c2,temp_c3;
	
	//Forward 0
	Xup=loclx_neighup[X][0];
	color_summ(temp_c0,in[Xup][0],in[Xup][2]);
	color_summ(temp_c1,in[Xup][1],in[Xup][3]);
	unsafe_su3_prod_color(out[X][0],conf[X][0],temp_c0);
	unsafe_su3_prod_color(out[X][1],conf[X][0],temp_c1);
	color_copy(out[X][2],out[X][0]);
	color_copy(out[X][3],out[X][1]);
      
	//Backward 0
	Xdw=loclx_neighdw[X][0];
	color_subt(temp_c0,in[Xdw][0],in[Xdw][2]);
	color_subt(temp_c1,in[Xdw][1],in[Xdw][3]);
	unsafe_su3_dag_prod_color(temp_c2,conf[Xdw][0],temp_c0);
	unsafe_su3_dag_prod_color(temp_c3,conf[Xdw][0],temp_c1);
	summassign_color(out[X][0],temp_c2);
	summassign_color(out[X][1],temp_c3);
	subtassign_color(out[X][2],temp_c2);
	subtassign_color(out[X][3],temp_c3);

	//Forward 1
	Xup=loclx_neighup[X][1];
	color_isumm(temp_c0,in[Xup][0],in[Xup][3]);
	color_isumm(temp_c1,in[Xup][1],in[Xup][2]);
	unsafe_su3_prod_color(temp_c2,conf[X][1],temp_c0);
	unsafe_su3_prod_color(temp_c3,conf[X][1],temp_c1);
	summassign_color(out[X][0],temp_c2);
	summassign_color(out[X][1],temp_c3);
	subtassign_icolor(out[X][2],temp_c3);
	subtassign_icolor(out[X][3],temp_c2);
	
	//Backward 1
	Xdw=loclx_neighdw[X][1];
	color_isubt(temp_c0,in[Xdw][0],in[Xdw][3]);
	color_isubt(temp_c1,in[Xdw][1],in[Xdw][2]);
	unsafe_su3_dag_prod_color(temp_c2,conf[Xdw][1],temp_c0);
	unsafe_su3_dag_prod_color(temp_c3,conf[Xdw][1],temp_c1);
	summassign_color(out[X][0],temp_c2);
	summassign_color(out[X][1],temp_c3);
	summassign_icolor(out[X][2],temp_c3);
	summassign_icolor(out[X][3],temp_c2);
    
	//Forward 2
	Xup=loclx_neighup[X][2];
	color_summ(temp_c0,in[Xup][0],in[Xup][3]);
	color_subt(temp_c1,in[Xup][1],in[Xup][2]);
	unsafe_su3_prod_color(temp_c2,conf[X][2],temp_c0);
	unsafe_su3_prod_color(temp_c3,conf[X][2],temp_c1);
	summassign_color(out[X][0],temp_c2);
	summassign_color(out[X][1],temp_c3);
	subtassign_color(out[X][2],temp_c3);
	summassign_color(out[X][3],temp_c2);
	
	//Backward 2
	Xdw=loclx_neighdw[X][2];
	color_subt(temp_c0,in[Xdw][0],in[Xdw][3]);
	color_summ(temp_c1,in[Xdw][1],in[Xdw][2]);
	unsafe_su3_dag_prod_color(temp_c2,conf[Xdw][2],temp_c0);
	unsafe_su3_dag_prod_color(temp_c3,conf[Xdw][2],temp_c1);
	summassign_color(out[X][0],temp_c2);
	summassign_color(out[X][1],temp_c3);
	summassign_color(out[X][2],temp_c3);
	subtassign_color(out[X][3],temp_c2);
    
	//Forward 3
	Xup=loclx_neighup[X][3];
	color_isumm(temp_c0,in[Xup][0],in[Xup][2]);
	color_isubt(temp_c1,in[Xup][1],in[Xup][3]);
	unsafe_su3_prod_color(temp_c2,conf[X][3],temp_c0);
	unsafe_su3_prod_color(temp_c3,conf[X][3],temp_c1);
	summassign_color(out[X][0],temp_c2);
	summassign_color(out[X][1],temp_c3);
	subtassign_icolor(out[X][2],temp_c2);
	summassign_icolor(out[X][3],temp_c3);
	
	//Backward 3
	Xdw=loclx_neighdw[X][3];
	color_isubt(temp_c0,in[Xdw][0],in[Xdw][2]);
	color_isumm(temp_c1,in[Xdw][1],in[Xdw][3]);
	unsafe_su3_dag_prod_color(temp_c2,conf[Xdw][3],temp_c0);
	unsafe_su3_dag_prod_color(temp_c3,conf[Xdw][3],temp_c1);
	summassign_color(out[X][0],temp_c2);
	summassign_color(out[X][1],temp_c3);
	summassign_icolor(out[X][2],temp_c2);
	subtassign_icolor(out[X][3],temp_c3);
	
	//for(int id=0;id<4;id++)
	//for(int ic=0;ic<3;ic++)
	//for(int ri=0;ri<2;ri++)
	//out[X][id][ic][ri]*=-0.5;
      }
}

void apply_Doe(spincolor *out,spincolor *in,quad_su3 *conf)
{apply_Doe_or_Deo(out,in,conf,0);}
void apply_Deo(spincolor *out,spincolor *in,quad_su3 *conf)
{apply_Doe_or_Deo(out,in,conf,1);}

void gamma5_ee_or_oo(spincolor *out,spincolor *in,int parity)
{
  for(int X=0;X<loc_vol;X++)
    if(loclx_parity[X]==parity)
      for(int ic=0;ic<3;ic++)
	for(int ri=0;ri<2;ri++)
	  {
	    out[X][0][ic][ri]=in[X][0][ic][ri];
	    out[X][1][ic][ri]=in[X][1][ic][ri];
	    out[X][2][ic][ri]=-in[X][2][ic][ri];
	    out[X][3][ic][ri]=-in[X][3][ic][ri];
	  }
}

void gamma5_ee(spincolor *out,spincolor *in)
{gamma5_ee_or_oo(out,in,0);}
void gamma5_oo(spincolor *out,spincolor *in)
{gamma5_ee_or_oo(out,in,1);}

void Qtm_minus_psi(spincolor *l,spincolor *k,quad_su3 *conf,spincolor *t0,spincolor *t1)
{
  apply_Deo(t1,k,conf);
  assign_mul_one_pm_imu_inv_oo(t1,t1,-1);
  apply_Doe(t0,t1,conf);
  assign_mul_one_pm_imu_sub_mul_gamma5_oo(l,k,t0,-1);
}

void Qtm_pm_psi(spincolor *l,spincolor *k,quad_su3 *conf,spincolor *t0,spincolor *t1,spincolor *t2)
{
  apply_Deo(t1,k,conf);
  assign_mul_one_pm_imu_inv_oo(t1,t1,-1);
  apply_Doe(t0,t1,conf);
  assign_mul_one_pm_imu_sub_mul_gamma5_oo(t0,k,t0,-1);
  
  apply_Deo(l,t0,conf);
  assign_mul_one_pm_imu_inv_oo(l,l,+1);
  apply_Doe(t1,l,conf);
  assign_mul_one_pm_imu_sub_mul_gamma5_oo(l,t0,t1,+1);
}

void cg_her(spincolor *sol,spincolor *source,spincolor *guess,quad_su3 *conf,int niter,int rniter,double residue)
{
  int riter=0;
  spincolor *p=allocate_spincolor(loc_vol+loc_bord,"p");
  spincolor *r=allocate_spincolor(loc_vol,"r");
  spincolor *s=allocate_spincolor(loc_vol,"s");
  spincolor *t1=allocate_spincolor(loc_vol+loc_bord,"t1"); //temporary for internal calculation of DD
  spincolor *t2=allocate_spincolor(loc_vol+loc_bord,"t2"); //temporary for internal calculation of DD
  spincolor *t3=allocate_spincolor(loc_vol+loc_bord,"t3");
  
  ///////////////// prepare the internal source /////////////////
  
  if(guess==NULL) memset(sol,0,sizeof(spincolor)*(loc_vol+loc_bord));
  else
    {
      communicate_lx_spincolor_borders(guess);
      memcpy(sol,guess,sizeof(spincolor)*(loc_vol+loc_bord));
    }

  //external loop, used if the internal exceed the maximal number of iterations
  double lambda; //(r_(k+1),r_(k+1))
  double source_norm;
  do
    {
      //calculate p0=r0=DD*sol_0 and delta_0=(p0,p0), performing global reduction and broadcast to all nodes
      double delta;
      {
	Qtm_pm_psi(s,sol,conf,t1,t2,t3);
	
        double loc_delta=0,loc_source_norm=0;
	for(int X=0;X<loc_vol;X++)
	  if(loclx_parity[X]==1)
	    for(int id=0;id<4;id++)
	      for(int ic=0;ic<3;ic++)
		for(int ri=0;ri<2;ri++)
		  {
		    double c1=source[X][id][ic][ri]-s[X][id][ic][ri];
		    p[X][id][ic][ri]=r[X][id][ic][ri]=c1;
		    if(riter==0) loc_source_norm+=source[X][id][ic][ri]*source[X][id][ic][ri];
		    loc_delta+=c1*c1;
		  }
        MPI_Allreduce(&loc_delta,&delta,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
        if(riter==0) MPI_Allreduce(&loc_source_norm,&source_norm,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
      }
      
      //main loop
      int iter=0;
      do
	{
	  double omega; //(r_k,r_k)/(p_k*DD*p_k)
	  
	  double alpha;
	  communicate_lx_spincolor_borders(p);
	  
	  //apply_Q2_RL(s,p,conf,t1,NULL,NULL,0);
	  Qtm_pm_psi(s,p,conf,t1,t2,t3);

	  double loc_alpha=0; //real part of the scalar product
	  for(int X=0;X<loc_vol;X++)
	    if(loclx_parity[X]==1)
	      for(int id=0;id<4;id++)
		for(int ic=0;ic<3;ic++)
		  for(int ri=0;ri<2;ri++)
		    loc_alpha+=s[X][id][ic][ri]*p[X][id][ic][ri];
	  MPI_Allreduce(&loc_alpha,&alpha,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
	  omega=delta/alpha;
	  printf("omega: %g\n",omega);

	  double loc_lambda=0;
	  for(int X=0;X<loc_vol;X++)
	    if(loclx_parity[X]==1)
	      for(int id=0;id<4;id++)
		for(int ic=0;ic<3;ic++)
		  for(int ri=0;ri<2;ri++)
		    {
		      if(debug3) printf("%d %d %d %d %g %g %g\n",X,id,ic,ri,s[X][id][ic][ri],r[X][id][ic][ri],p[X][id][ic][ri]);
		      sol[X][id][ic][ri]+=omega*p[X][id][ic][ri];    //sol_(k+1)=x_k+omega*p_k
		      double c1=r[X][id][ic][ri]-omega*s[X][id][ic][ri];//r_(k+1)=x_k-omega*pk
		      r[X][id][ic][ri]=c1;
		      loc_lambda+=c1*c1;
		    }
	  MPI_Allreduce(&loc_lambda,&lambda,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
	  
	  double gammag=lambda/delta; //(r_(k+1),r_(k+1))/(r_k,r_k)
	  delta=lambda;
	  
	  //p_(k+1)=r_(k+1)+gammag*p_k
	  for(int X=0;X<loc_vol;X++)
	    if(loclx_parity[X]==1)
	      for(int id=0;id<4;id++)
		for(int ic=0;ic<3;ic++)
		  for(int ri=0;ri<2;ri++)
		    {
		      if(debug3) printf("r=%g p=%g\n",r[X][id][ic][ri],p[X][id][ic][ri]);
		      p[X][id][ic][ri]=r[X][id][ic][ri]+gammag*p[X][id][ic][ri];
		    }
	  
	  iter++;
	  
	  if(rank==0 && debug) printf("iter %d residue %g\n",iter,lambda);
	}
      while(lambda>(residue*source_norm) && iter<niter);
      
      //last calculation of residual, in the case iter>niter
      communicate_lx_spincolor_borders(sol);

      Qtm_pm_psi(s,sol,conf,t1,t2,t3);
      {
        double loc_lambda=0;
	for(int X=0;X<loc_vol;X++)
	  if(loclx_parity[X]==1)
	    for(int id=0;id<4;id++)
	      for(int ic=0;ic<3;ic++)
		for(int ri=0;ri<2;ri++)
		  {
		    double c1=source[X][id][ic][ri]-s[X][id][ic][ri];
		    loc_lambda+=c1*c1;
		  }
	if(rank_tot>0) MPI_Allreduce(&loc_lambda,&lambda,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
	else lambda=loc_lambda;
      }

      riter++;
    }
  while(lambda>(residue*source_norm) && riter<rniter);
  
  free(s);
  free(p);
  free(r);
  free(t1);
  free(t2);
  free(t3);
}

void inv_Q_eo_impr_cg(spincolor *Even_new,spincolor *Odd_new,spincolor *Even,spincolor *Odd,spincolor *guess,quad_su3 *conf,int niter,int rniter,double residue)
 {
  spincolor *t0=allocate_spincolor(loc_vol+loc_bord,"t0");
  spincolor *t1=allocate_spincolor(loc_vol+loc_bord,"t1");

  //create the source for the OO
  assign_mul_one_pm_imu_inv_ee(Even_new,Even,1);
  apply_Doe(t0,Even_new,conf);
  assign_mul_add_r_oo(t0,1,Odd);
  gamma5_oo(t0,t0);

  cg_her(Odd_new,t0,guess,conf,niter,rniter,residue);
  
  //now remove the additional Q
  Qtm_minus_psi(Odd_new,Odd_new,conf,t0,t1);
  
  //reconstruct even sites
  apply_Deo(t0,Odd_new,conf);
  assign_mul_one_pm_imu_inv_ee(t0,t0,1);
  assign_add_mul_r_ee(Even_new,t0,1);
  
  free(t0);
  free(t1);
}

void mega_wrap(spincolor *sol,spincolor *source,spincolor *guess,quad_su3 *conf,int niter,int rniter,double residue)
{
  spincolor *Even=allocate_spincolor(loc_vol,"Even");
  spincolor *Odd=allocate_spincolor(loc_vol,"Odd");
  spincolor *Even_new=allocate_spincolor(loc_vol,"Even_new");
  spincolor *Odd_new=allocate_spincolor(loc_vol,"Odd_new");
  
  //split in even and odd
  for(int X=0;X<loc_vol;X++)
    if(loclx_parity[X]==0) memcpy(Even[X],source[X],sizeof(spincolor));
    else                   memcpy(Odd[X],source[X],sizeof(spincolor));
  
  inv_Q_eo_impr_cg(Even_new,Odd_new,Even,Odd,guess,conf,niter,rniter,residue);
  
  //repaste even and odd
  for(int X=0;X<loc_vol;X++)
    if(loclx_parity[X]==0) memcpy(sol[X],Even_new[X],sizeof(spincolor));
    else                   memcpy(sol[X],Odd_new[X],sizeof(spincolor));
  
  free(Odd_new);
  free(Even_new);
  free(Odd);
  free(Even);
}

int main(int narg,char **arg)
{
  //basic mpi initialization
  init_appretto();

  if(narg<2 && rank==0)
    {
      fprintf(stderr,"Use: %s input_file\n",arg[0]);
      fflush(stderr);
      MPI_Abort(MPI_COMM_WORLD,1);
    }

  open_input(arg[1]);

  read_str_int("L",&(glb_size[1]));
  read_str_int("T",&(glb_size[0]));

  read_str_double("m",&(mu));
  read_str_double("kappa",&(kappa));
  muh=2*kappa*mu;
  
  //Init the MPI grid 
  init_grid();
  set_eo_geometry();
  
  //Initialize the gauge configuration and read the path
  quad_su3 *conf=(quad_su3*)malloc(sizeof(quad_su3)*(loc_vol+loc_bord));
  char gauge_file[1024];
  read_str_str("GaugeConf",gauge_file,1024);
  
  //read the thetas in multiples of 2*pi
  double theta[4];
  read_str_double("ThetaTXYZ",&(theta[0]));
  read_double(&(theta[1]));
  read_double(&(theta[2]));
  read_double(&(theta[3]));
  if(rank==0)
    printf("Thetas: %f %f %f %f\n",theta[0],theta[1],theta[2],theta[3]);

  //load the configuration, put boundaries condition and communicate borders
  read_local_gauge_conf(conf,gauge_file);
  put_boundaries_conditions(conf,theta,1,0);
  communicate_gauge_borders(conf);

  //initialize source to delta
  spincolor *source=allocate_spincolor(loc_vol+loc_bord,"source");
  memset(source,0,sizeof(spincolor)*loc_vol);
  if(rank==0) source[0][0][0][0]=1;

  double residue;
  read_str_double("Residue",&residue);
  int nitermax;
  read_str_int("NiterMax",&nitermax);
  
  communicate_lx_spincolor_borders(source);

  //initialize solution
  spincolor *solutionQ2=(spincolor*)malloc(sizeof(spincolor)*(loc_vol+loc_bord));
  spincolor *solutionQ=(spincolor*)malloc(sizeof(spincolor)*(loc_vol+loc_bord));

  close_input();
  
  ///////////////////////////////////////////
  
  spincolor *temp0=allocate_spincolor(loc_vol+loc_bord,"temp0");
  spincolor *tempe=allocate_spincolor(loc_vol+loc_bord,"tempe");
  spincolor *tempo=allocate_spincolor(loc_vol+loc_bord,"tempe");
  spincolor *temp1=allocate_spincolor(loc_vol+loc_bord,"temp1");
  //spincolor *temp2=allocate_spincolor(loc_vol+loc_bord,"temp2");

  for(int X=0;X<loc_vol;X++)
    for(int id=0;id<4;id++)
      for(int ic=0;ic<3;ic++)
	for(int ri=0;ri<2;ri++)
	  temp0[X][id][ic][ri]=(double)rand()/RAND_MAX;
  
  apply_Deo(tempe,temp0,conf);
  apply_Doe(tempo,temp0,conf);

  for(int X=0;X<loc_vol;X++)
    for(int id=0;id<4;id++)
      for(int ic=0;ic<3;ic++)
	for(int ri=0;ri<2;ri++)
	  {
	    if(loclx_parity[X]==0) temp1[X][id][ic][ri]=tempe[X][id][ic][ri];
	    else                   temp1[X][id][ic][ri]=tempo[X][id][ic][ri];
	    
	    //printf("%d %d %d %d %g %g\n",X,id,ic,ri,temp1[X][id][ic][ri],temp2[X][id][ic][ri]);
	  }

  //take initial time                                                                                                        
  double tic;
  MPI_Barrier(cart_comm);
  tic=MPI_Wtime();
  inv_Q2_cg(solutionQ,source,NULL,conf,kappa,mu,nitermax,1,residue);
  apply_Q(solutionQ2,solutionQ,conf,kappa,mu);
  mega_wrap(solutionQ,source,NULL,conf,nitermax,1,residue);

  MPI_Barrier(cart_comm);
  double tac=MPI_Wtime();
  if(rank==0)
    printf("\nTotal time elapsed: %f s\n",tac-tic);

  spincolor *source_reco=allocate_spincolor(loc_vol,"source_reco");
  
  apply_Q(source_reco,solutionQ,conf,kappa,mu);
  
  //printing
  double truered,loc_truered=0;
  for(int loc_site=0;loc_site<loc_vol;loc_site++)
    for(int id=0;id<4;id++)
      for(int ic=0;ic<3;ic++)
	{
	  for(int ri=0;ri<2;ri++)
	    printf("%g %g\n",solutionQ2[loc_site][id][ic][ri],solutionQ[loc_site][id][ic][ri]);
	  double tempr=source[loc_site][id][ic][0]-source_reco[loc_site][id][ic][0];
	  double tempi=source[loc_site][id][ic][1]-source_reco[loc_site][id][ic][1];
	  loc_truered+=tempr*tempr+tempi*tempi;
	}

  if(rank_tot>0) MPI_Allreduce(&loc_truered,&truered,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);                                 
  else truered=loc_truered;

  if(rank==0) printf("Residue: %g\n",truered); 

  ///////////////////////////////////////////

  /*
  theta[0]=0;
  theta[1]=theta[2]=theta[3]=0.1;
  put_boundaries_conditions(conf,theta,1,0);
  communicate_gauge_borders(conf);
  inv_Q2_cg(solution,source,solution,conf,1000,1,1.e-6);
  */

  free(solutionQ);
  free(solutionQ2);
  free(source);
  free(source_reco);

  free(conf);

  ///////////////////////////////////////////

  close_appretto();

  return 0;
}
