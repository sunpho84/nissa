#include "nissa.h"

const int debug_lvl3=0;

spincolor *glb_source;
quad_su3 *conf;
int nitermax;
double residue;

void gamma5(spincolor *out,spincolor *in)
{
  for(int ivol=0;ivol<loc_vol;ivol++)
    for(int ic=0;ic<3;ic++)
      for(int ri=0;ri<2;ri++)
	{
	  for(int id=0;id<2;id++) out[ivol][id][ic][ri]=+in[ivol][id][ic][ri];
	  for(int id=2;id<4;id++) out[ivol][id][ic][ri]=-in[ivol][id][ic][ri];
	}
}

//apply D_eo or oe
void Doe_or_Deo(spincolor *out,spincolor *in,quad_su3 *conf,int parity)
{
  int Xup,Xdw;

  if(debug_lvl3) printf(" m2Doe_or_m2Deo parity: %d\n",parity);
  for(int X=0;X<loc_vol;X++)
    if(loclx_parity[X]!=parity) //loop is on sink
      {
	memset(out[X],0,sizeof(spincolor));
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
        color_summassign(out[X][0],temp_c2);
        color_summassign(out[X][1],temp_c3);
        color_subtassign(out[X][2],temp_c2);
        color_subtassign(out[X][3],temp_c3);

        //Forward 1
        Xup=loclx_neighup[X][1];
        color_isumm(temp_c0,in[Xup][0],in[Xup][3]);
        color_isumm(temp_c1,in[Xup][1],in[Xup][2]);
        unsafe_su3_prod_color(temp_c2,conf[X][1],temp_c0);
        unsafe_su3_prod_color(temp_c3,conf[X][1],temp_c1);
        color_summassign(out[X][0],temp_c2);
        color_summassign(out[X][1],temp_c3);
        color_isubtassign(out[X][2],temp_c3);
        color_isubtassign(out[X][3],temp_c2);
        
        //Backward 1
        Xdw=loclx_neighdw[X][1];
        color_isubt(temp_c0,in[Xdw][0],in[Xdw][3]);
        color_isubt(temp_c1,in[Xdw][1],in[Xdw][2]);
        unsafe_su3_dag_prod_color(temp_c2,conf[Xdw][1],temp_c0);
        unsafe_su3_dag_prod_color(temp_c3,conf[Xdw][1],temp_c1);
        color_summassign(out[X][0],temp_c2);
        color_summassign(out[X][1],temp_c3);
        color_isummassign(out[X][2],temp_c3);
        color_isummassign(out[X][3],temp_c2);
    
        //Forward 2
        Xup=loclx_neighup[X][2];
        color_summ(temp_c0,in[Xup][0],in[Xup][3]);
        color_subt(temp_c1,in[Xup][1],in[Xup][2]);
        unsafe_su3_prod_color(temp_c2,conf[X][2],temp_c0);
        unsafe_su3_prod_color(temp_c3,conf[X][2],temp_c1);
        color_summassign(out[X][0],temp_c2);
        color_summassign(out[X][1],temp_c3);
        color_subtassign(out[X][2],temp_c3);
        color_summassign(out[X][3],temp_c2);
        
        //Backward 2
        Xdw=loclx_neighdw[X][2];
        color_subt(temp_c0,in[Xdw][0],in[Xdw][3]);
        color_summ(temp_c1,in[Xdw][1],in[Xdw][2]);
        unsafe_su3_dag_prod_color(temp_c2,conf[Xdw][2],temp_c0);
        unsafe_su3_dag_prod_color(temp_c3,conf[Xdw][2],temp_c1);
        color_summassign(out[X][0],temp_c2);
        color_summassign(out[X][1],temp_c3);
        color_summassign(out[X][2],temp_c3);
        color_subtassign(out[X][3],temp_c2);
    
        //Forward 3
        Xup=loclx_neighup[X][3];
        color_isumm(temp_c0,in[Xup][0],in[Xup][2]);
        color_isubt(temp_c1,in[Xup][1],in[Xup][3]);
        unsafe_su3_prod_color(temp_c2,conf[X][3],temp_c0);
        unsafe_su3_prod_color(temp_c3,conf[X][3],temp_c1);
        color_summassign(out[X][0],temp_c2);
        color_summassign(out[X][1],temp_c3);
        color_isubtassign(out[X][2],temp_c2);
        color_isummassign(out[X][3],temp_c3);
        
        //Backward 3
        Xdw=loclx_neighdw[X][3];
        color_isubt(temp_c0,in[Xdw][0],in[Xdw][2]);
        color_isumm(temp_c1,in[Xdw][1],in[Xdw][3]);
        unsafe_su3_dag_prod_color(temp_c2,conf[Xdw][3],temp_c0);
        unsafe_su3_dag_prod_color(temp_c3,conf[Xdw][3],temp_c1);
        color_summassign(out[X][0],temp_c2);
        color_summassign(out[X][1],temp_c3);
        color_isummassign(out[X][2],temp_c2);
        color_isubtassign(out[X][3],temp_c3);

	for(int id=0;id<4;id++)
          for(int ic=0;ic<3;ic++)
            for(int ri=0;ri<2;ri++)
		out[X][id][ic][ri]*=-0.5;
      }
}

void direct_invert(spincolor *solution,spincolor *source,double mu,double kappa)
{
  spincolor *temp_source=nissa_malloc("temp_source",loc_vol,spincolor);
  spincolor *solutionQ2=nissa_malloc("solutionQ2",loc_vol,spincolor);

  gamma5(temp_source,source);
  inv_tmQ2_cg(solutionQ2,temp_source,NULL,conf,kappa,mu,nitermax,1,residue);
  apply_tmQ(solution,solutionQ2,conf,kappa,-mu);

  nissa_free(temp_source);
  nissa_free(solutionQ2);
}

void dir_inv_Dee_Doo(spincolor *out,spincolor *in,double mu,double kappa,int parity,int dir_inv)
{
  double a=1/(2*kappa);
  double b=mu;
  complex z;
  if(dir_inv==0)
    {
      z[0]=a;
      z[1]=b;
    }
  else
    {
      double nrm=a*a+b*b;
      z[0]=+a/nrm;
      z[1]=-b/nrm;
    }
  
  spincolor temp;
  for(int X=0;X<loc_vol;X++)
    if(loclx_parity[X]==parity)
      for(int ic=0;ic<3;ic++)
	{
	  for(int id=0;id<2;id++) unsafe_complex_prod(temp[id][ic],in[X][id][ic],z);
	  for(int id=2;id<4;id++) unsafe_complex_conj2_prod(temp[id][ic],in[X][id][ic],z);
	  memcpy(out[X],temp,sizeof(spincolor));
	}
}

void dir_Dee_Doo(spincolor *out,spincolor *in,double mu,double kappa,int parity)
{dir_inv_Dee_Doo(out,in,mu,kappa,parity,0);}
void inv_Dee_Doo(spincolor *out,spincolor *in,double mu,double kappa,int parity)
{dir_inv_Dee_Doo(out,in,mu,kappa,parity,1);}

void inv_Dee(spincolor *out,spincolor *in,double mu,double kappa){inv_Dee_Doo(out,in,mu,kappa,0);}
void inv_Doo(spincolor *out,spincolor *in,double mu,double kappa){inv_Dee_Doo(out,in,mu,kappa,1);}

void Dee(spincolor *out,spincolor *in,double mu,double kappa)    {dir_Dee_Doo(out,in,mu,kappa,0);}
void Doo(spincolor *out,spincolor *in,double mu,double kappa)    {dir_Dee_Doo(out,in,mu,kappa,1);}
void Doe(spincolor *out,spincolor *in,quad_su3*conf){Doe_or_Deo(out,in,conf,0);}
void Deo(spincolor *out,spincolor *in,quad_su3*conf){Doe_or_Deo(out,in,conf,1);}

void K(spincolor *out,spincolor *in,quad_su3* conf,double mu,double kappa)
{
  spincolor *temp=nissa_malloc("temp",loc_vol,spincolor);

  Deo(out,in,conf);
  inv_Dee(temp,out,mu,kappa);
  Doe(out,temp,conf);  
  inv_Doo(temp,out,mu,kappa);
  Doo(temp,in,mu,kappa);

  for(int ivol=0;ivol<loc_vol;ivol++)
    if(loclx_parity[ivol]==1)
      for(int id=0;id<4;id++)
	for(int ic=0;ic<3;ic++)
	  for(int ri=0;ri<2;ri++)
	    out[ivol][id][ic][ri]=temp[ivol][id][ic][ri]-out[ivol][id][ic][ri];

  gamma5(out,out);
  
  nissa_free(temp);
}

void K2(spincolor *out,spincolor *in,quad_su3 *conf,double mu,double kappa)
{
  memset(out,0,sizeof(spincolor)*loc_vol);
  
  spincolor *temp=nissa_malloc("temp",loc_vol,spincolor);
  K(temp,in,conf,-mu,kappa);
  K(out,temp,conf,+mu,kappa);
  nissa_free(temp);
}

void inv_K(spincolor *sol,spincolor *source,double mu,double kappa)
{
  int niter=nitermax;
  int riter=0;
  int rniter=5;
  spincolor *guess=NULL;
  spincolor *p=nissa_malloc("p",loc_vol+bord_vol,spincolor);
  spincolor *r=nissa_malloc("r",loc_vol,spincolor);
  spincolor *s=nissa_malloc("s",loc_vol,spincolor);
  
  ///////////////// prepare the internal source /////////////////
  
  if(guess==NULL) memset(sol,0,sizeof(spincolor)*(loc_vol+bord_vol));
  else
    {
      communicate_lx_spincolor_borders(guess);
      memcpy(sol,guess,sizeof(spincolor)*(loc_vol+bord_vol));
    }

  //external loop, used if the internal exceed the maximal number of iterations
  double lambda; //(r_(k+1),r_(k+1))
  double source_norm;
  do
    {
      //calculate p0=r0=DD*sol_0 and delta_0=(p0,p0), performing global reduction and broadcast to all nodes
      double delta;
      {
	K2(s,sol,conf,mu,kappa);
	
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
	  
	  K2(s,p,conf,mu,kappa);
	  
	  double loc_alpha=0; //real part of the scalar product
	  for(int X=0;X<loc_vol;X++)
	    if(loclx_parity[X]==1)
	      for(int id=0;id<4;id++)
		for(int ic=0;ic<3;ic++)
		  for(int ri=0;ri<2;ri++)
		    loc_alpha+=s[X][id][ic][ri]*p[X][id][ic][ri];
	  MPI_Allreduce(&loc_alpha,&alpha,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
	  omega=delta/alpha;
	  //printf("omega: %g\n",omega);
	  
	  double loc_lambda=0;
	  for(int X=0;X<loc_vol;X++)
	    if(loclx_parity[X]==1)
	      for(int id=0;id<4;id++)
		for(int ic=0;ic<3;ic++)
		  for(int ri=0;ri<2;ri++)
		    {
		      if(debug_lvl3) printf("%d %d %d %d %g %g %g\n",X,id,ic,ri,s[X][id][ic][ri],r[X][id][ic][ri],p[X][id][ic][ri]);
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
		      if(debug_lvl3) printf("r=%g p=%g\n",r[X][id][ic][ri],p[X][id][ic][ri]);
		      p[X][id][ic][ri]=r[X][id][ic][ri]+gammag*p[X][id][ic][ri];
		    }
	  
	  iter++;
	  
	  if(rank==0 && debug_lvl) printf("eo iter %d residue %g\n",iter,lambda);
	}
      while(lambda>(residue*source_norm) && iter<niter);
      
      //last calculation of residual, in the case iter>niter
      communicate_lx_spincolor_borders(sol);

      K2(s,sol,conf,mu,kappa);
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
	if(nissa_nranks>0) MPI_Allreduce(&loc_lambda,&lambda,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
	else lambda=loc_lambda;
      }

      riter++;
    }
  while(lambda>(residue*source_norm) && riter<rniter);

  K(s,sol,conf,-mu,kappa);
  memcpy(sol,s,sizeof(spincolor)*loc_vol);
  
  nissa_free(s);
  nissa_free(p);
  nissa_free(r);
}

void improved_invert(spincolor *solution,spincolor *chi,double mu,double kappa)
{
  spincolor *varphi=nissa_malloc("varphi",loc_vol,spincolor);

  spincolor *temp=nissa_malloc("temp",loc_vol,spincolor);
  inv_Dee(temp,chi,mu,kappa);
  Doe(varphi,temp,conf);
  nissa_free(temp);
  
  for(int ivol=0;ivol<loc_vol;ivol++)
    if(loclx_parity[ivol]==1)
      for(int id=0;id<2;id++)
	for(int ic=0;ic<3;ic++)
	  for(int ri=0;ri<2;ri++)
	    {    
	      varphi[ivol][id  ][ic][ri]=+chi[ivol][id  ][ic][ri]-varphi[ivol][id  ][ic][ri];
	      varphi[ivol][id+2][ic][ri]=-chi[ivol][id+2][ic][ri]+varphi[ivol][id+2][ic][ri];
	    }
  
  inv_K(solution,varphi,mu,kappa);

  Deo(varphi,solution,conf);
  for(int ivol=0;ivol<loc_vol;ivol++)
    if(loclx_parity[ivol]==0)
      for(int id=0;id<4;id++)
	for(int ic=0;ic<3;ic++)
	  for(int ri=0;ri<2;ri++)
	    varphi[ivol][id][ic][ri]=chi[ivol][id][ic][ri]-varphi[ivol][id][ic][ri];
  inv_Dee(solution,varphi,mu,kappa);

  nissa_free(varphi);
}

void init(char *input_path,double *mu,double *kappa)
{
  open_input(input_path);

  int L,T;
  read_str_int("L",&L);
  read_str_int("T",&T);

  read_str_double("m",mu);
  read_str_double("kappa",kappa);
  
  //Init the MPI grid 
  init_grid(T,L);
  //set_eo_geometry();
  
  //Initialize the gauge configuration and read the path
  conf=nissa_malloc("conf",loc_vol+bord_vol,quad_su3);
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
  read_ildg_gauge_conf(conf,gauge_file);
  put_boundaries_conditions(conf,theta,1,0);
  communicate_lx_quad_su3_borders(conf);

  //initialize source to delta
  glb_source=nissa_malloc("source",loc_vol+bord_vol,spincolor);
  memset(glb_source,0,sizeof(spincolor)*loc_vol);
  if(rank==0)
    {
      glb_source[1][0][0][0]=1;
      //glb_source[0][0][0][0]=1;
    }
  communicate_lx_spincolor_borders(glb_source);
  
  read_str_double("Residue",&residue);
  read_str_int("NiterMax",&nitermax);

  close_input();
}

int main(int narg,char **arg)
{
  double mu;
  double kappa;
  
  //basic mpi initialization
  init_nissa();

  if(narg<2 && rank==0)
    {
      fprintf(stderr,"Use: %s input_file\n",arg[0]);
      fflush(stderr);
      MPI_Abort(MPI_COMM_WORLD,1);
    }

  init(arg[1],&mu,&kappa);

  //initialize solution
  spincolor *solution_direct=nissa_malloc("direct_sol",loc_vol+bord_vol,spincolor);
  spincolor *solution_improved=nissa_malloc("direct_sol",loc_vol+bord_vol,spincolor);
  
  ///////////////////////////////////////////
  
  direct_invert(solution_direct,glb_source,mu,kappa);
  improved_invert(solution_improved,glb_source,mu,kappa);

  for(int ivol=0;ivol<loc_vol;ivol++)
    for(int id=0;id<4;id++)
      for(int ic=0;ic<3;ic++)
	for(int ri=0;ri<2;ri++)
	  printf("%lg %lg\n",solution_direct[ivol][id][ic][ri],solution_improved[ivol][id][ic][ri]);
  
  ///////////////////////////////////////////

  nissa_free(solution_direct);
  nissa_free(solution_improved);
  nissa_free(glb_source);
  
  nissa_free(conf);

  ///////////////////////////////////////////

  close_nissa();

  return 0;
}
