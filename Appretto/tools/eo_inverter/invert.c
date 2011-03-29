#include <mpi.h>
#include <lemon.h>

#include "appretto.h"

const int debug3=0;

void inv_Doo_or_Dee(spincolor *psi,spincolor *source,double kappa,double mu,int parity)
{
  double muh=2*kappa*mu;
  complex factp={1/(1+muh*muh),muh/(1+muh*muh)};
  complex factm={factp[0],-factm[1]};
  
  for(int X=0;X<loc_vol;X++)
    if(loclx_parity[X]==parity)
      {
	for(int id=0;id<2;id++)
	  for(int ic=0;ic<3;ic++)
	    {
	      unsafe_complex_prod(psi[X][id][ic],source[X][id][ic],factm);
	      unsafe_complex_prod(psi[X][id+2][ic],source[X][id+2][ic],factp);
	    }
      }
}

void inv_Doo(spincolor *psi,spincolor *source,double kappa,double mu)
{
  inv_Doo_or_Dee(psi,source,kappa,mu,1);
}

void inv_Dee(spincolor *psi,spincolor *source,double kappa,double mu)
{
  inv_Doo_or_Dee(psi,source,kappa,mu,0);
}

//apply D_eo or oe
void apply_Doe_or_Deo(spincolor *out,spincolor *in,quad_su3 *conf,double kappa,int parity)
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
	
	for(int id=0;id<4;id++)
	  for(int ic=0;ic<3;ic++)
	    for(int ri=0;ri<2;ri++)
	      out[X][id][ic][ri]*=-kappa;
      }
}

/*
void apply_Doe_or_Deo(spincolor *out,spincolor *in,quad_su3 *conf,int parity)
{
  apply_m2Doe_or_m2Deo(out,in,conf,parity);
  if(debug3) printf(" Doe_or_doe parity: %d\n",parity);
  for(int X=0;X<loc_vol;X++)
    if(loclx_parity[X]!=parity) //this loop is on sink
      for(int id=0;id<4;id++)
	for(int ic=0;ic<3;ic++)
	  for(int ri=0;ri<2;ri++)
	    {
	      if(debug3) printf("X=%d\n",X);
	      out[X][id][ic][ri]*=-0.5;
	    }
}

void apply_Doe(spincolor *out,spincolor *in,quad_su3 *conf)
{
  apply_Doe_or_Deo(out,in,conf,0);
}

void apply_Deo(spincolor *out,spincolor *in,quad_su3 *conf)
{
  apply_Doe_or_Deo(out,in,conf,1);
}
*/

void apply_Doe(spincolor *out,spincolor *in,double kappa,quad_su3 *conf)
{
  apply_Doe_or_Deo(out,in,conf,kappa,0);
}

void apply_Deo(spincolor *out,spincolor *in,double kappa,quad_su3 *conf)
{
  apply_Doe_or_Deo(out,in,conf,kappa,1);
}

void apply_Q2_oo_impr(spincolor *s,spincolor *p,quad_su3 *conf,double kappa,double m,spincolor *t1,spincolor *t2)
{
  apply_Deo(t1,p,kappa,conf);
  inv_Dee(t2,t1,kappa,m);
  communicate_lx_spincolor_borders(t2);
  apply_Doe(t1,t2,kappa,conf);
  inv_Doo(t2,t1,kappa,m);
  
  for(int X=0;X<loc_vol;X++)
    if(loclx_parity[X]==1)
      for(int id=0;id<4;id++)
	for(int ic=0;ic<3;ic++)
	  for(int ri=0;ri<2;ri++)
	    s[X][id][ic][ri]=p[X][id][ic][ri]-t2[X][id][ic][ri];
}

void inv_Q_eo_impr_cg(spincolor *sol,spincolor *or_source,spincolor *guess,quad_su3 *conf,double kappa,double m,int niter,int rniter,double residue)
{
  int riter=0;
  spincolor *source=allocate_spincolor(loc_vol+loc_bord,"source");
  spincolor *s=allocate_spincolor(loc_vol,"s");
  spincolor *p=allocate_spincolor(loc_vol+loc_bord,"p");
  spincolor *r=allocate_spincolor(loc_vol,"r");
  spincolor *t1=allocate_spincolor(loc_vol+loc_bord,"t1"); //temporary for internal calculation of DD
  spincolor *t2=allocate_spincolor(loc_vol+loc_bord,"t2"); //temporary for internal calculation of DD

  ///////////////// prepare the internal source /////////////////
  communicate_lx_spincolor_borders(or_source);
  //even part
  inv_Dee(source,or_source,kappa,m);
  //odd part with the 20% claimed improvement
  communicate_lx_spincolor_borders(source);
  apply_Doe(t1,source,kappa,conf);

  for(int X=0;X<loc_vol;X++)
    if(loclx_parity[X]==1)
      for(int id=0;id<4;id++)
	for(int ic=0;ic<3;ic++)
	  for(int ri=0;ri<2;ri++)
	    p[X][id][ic][ri]=or_source[X][id][ic][ri]+t1[X][id][ic][ri];
  inv_Doo(source,p,kappa,m); //20% claimed improvement

  communicate_lx_spincolor_borders(source);
  
  //reset the solution
  memset(sol,0,sizeof(spincolor)*(loc_vol+loc_bord));
  
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
	apply_Q2_oo_impr(s,sol,conf,kappa,m,t1,t2);

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
	  
	  apply_Q2_oo_impr(s,p,conf,kappa,m,t1,t2);

	  double loc_alpha=0; //real part of the scalar product
	  for(int X=0;X<loc_vol;X++)
	    if(loclx_parity[X]==1)
	      for(int id=0;id<4;id++)
		for(int ic=0;ic<3;ic++)
		  for(int ri=0;ri<2;ri++)
		    loc_alpha+=s[X][id][ic][ri]*p[X][id][ic][ri];
	  if(rank_tot>0) MPI_Allreduce(&loc_alpha,&alpha,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
	  else alpha=loc_alpha;
	  omega=delta/alpha;
	  
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
	  if(rank_tot>0) MPI_Allreduce(&loc_lambda,&lambda,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
	  else lambda=loc_lambda;
	  
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

      apply_Q2_oo_impr(s,sol,conf,kappa,m,t1,t2);
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
  
  communicate_lx_spincolor_borders(sol);
  apply_Deo(t1,sol,kappa,conf);
  inv_Dee(p,t1,kappa,m);
  for(int X=0;X<loc_vol;X++)
    if(loclx_parity[X]==0)
      for(int id=0;id<4;id++)
	for(int ic=0;ic<3;ic++)
	  for(int ri=0;ri<2;ri++)
	    sol[X][id][ic][ri]=source[X][id][ic][ri]+0.5*p[X][id][ic][ri];
  
  free(s);
  free(p);
  free(r);
  free(t1);
  free(t2);
  free(source);
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

  double m;
  double kappa;
  read_str_double("m",&(m));
  read_str_double("kappa",&(kappa));

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
  spincolor *source=(spincolor*)malloc(sizeof(spincolor)*(loc_vol+loc_bord));
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
  
  spincolor *tempe=allocate_spincolor(loc_vol+loc_bord,"tempe");
  spincolor *tempo=allocate_spincolor(loc_vol+loc_bord,"tempe");
  spincolor *temp1=allocate_spincolor(loc_vol+loc_bord,"temp");
  
  apply_Deo(tempe,source,kappa,conf);
  apply_Doe(tempo,source,kappa,conf);

  apply_Deo(temp1,source,kappa,conf);
  apply_Doe(temp1,source,kappa,conf);
  
  for(int X=0;X<loc_vol;X++)
    for(int id=0;id<4;id++)
      for(int ic=0;ic<3;ic++)
	for(int ri=0;ri<2;ri++)
	  {
	    if(loclx_parity[X]==0) temp1[X][ic][id][ri]-=tempe[X][id][ic][ri];
	    else                   temp1[X][ic][id][ri]-=tempo[X][id][ic][ri];
	    
	    if(temp1[X][id][ic][ri]) printf("par=%d X=%d id=%d ic=%d ri=%d %g\n",loclx_parity[X],X,id,ic,ri,temp1[X][id][ic][ri]);
	  }

  //take initial time                                                                                                        
  double tic;
  MPI_Barrier(cart_comm);
  tic=MPI_Wtime();
  inv_Q_eo_impr_cg(solutionQ,source,NULL,conf,kappa,m,nitermax,1,residue);
  //inv_Q2_cg(solutionQ2,source,NULL,conf,kappa,m,nitermax,1,residue);

  MPI_Barrier(cart_comm);
  double tac=MPI_Wtime();
  if(rank==0)
    printf("\nTotal time elapsed: %f s\n",tac-tic);

  spincolor *source_reco=(spincolor*)malloc(sizeof(spincolor)*loc_vol);
  
  apply_Q2(source_reco,solutionQ,conf,kappa,m,NULL,NULL,NULL);
  
  //printing
  double truered,loc_truered=0;
  for(int loc_site=0;loc_site<loc_vol;loc_site++)
    for(int id=0;id<4;id++)
      for(int ic=0;ic<3;ic++)
	{
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
  inv_Q2_cg(solution,source,solution,conf,kappa,m,1000,1,1.e-6);
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
