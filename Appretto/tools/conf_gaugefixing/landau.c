#include "appretto.h"

double alpha;
double precision;

//return the log of the factorial of n
double lfact(double n)
{
  double log_factn=0;
  for(int i=1;i<=n;i++) log_factn+=logf(i);
  return log_factn;
}

//overrelax the transformation
void overrelax(su3 *out,su3 *in,int N,double omega)
{
  //find coefficients
  double coef[N+1];
  for(int n=0;n<=N;n++) coef[n]=exp(lgamma(omega+1)-lgamma(omega+1-n)-lfact(n));
  
  //loop over volume
  for(int ivol=0;ivol<loc_vol;ivol++)
    {
      su3 t,o,w;

      su3_summ_real(w,in[ivol],-1); //subtract 1 from in[ivol]
      su3_copy(t,w);                //init t
      su3_put_to_id(o);             //output init
      
      //compute various powers
      for(int n=1;n<=N;n++)
	{ //summ 
	  for(int i=0;i<18;i++) ((double*)o)[i]+=coef[n]*((double*)t)[i];
	  safe_su3_prod_su3(t,t,w); //next power of w
	}
      
      //unitarize
      su3_unitarize(out[ivol],o);
    }
}

//apply the fast fourier acceleration
void fast_fourier_accelerate_fixing(su3 *g)
{  
  fft4d((complex*)g,(complex*)g,9,1);
  
  //apply the conditioning
  for(int ivol=0;ivol<loc_vol;ivol++)
    {
      //compute p2
      double p2=0;
      for(int mu=0;mu<4;mu++)
	{
	  int ix=glb_coord_of_loclx[ivol][mu];
	  double pmu=M_PI*ix/glb_size[mu];
	  double sinpmu=sin(pmu);
	  p2+=sinpmu*sinpmu;
	}
      p2/=4;
      
      //apply
      if(p2!=0)
	for(int i=0;i<3;i++)
	  for(int j=0;j<3;j++)
	    for(int ri=0;ri<2;ri++)
	      g[ivol][i][j][ri]/=p2;
    }
  
    fft4d((complex*)g,(complex*)g,9,-1);
}

//compute the steepest descent gauge fixing transformation
void find_steepest_descent_fixing(su3 *g,quad_su3 *conf,double alpha)
{
  memset(g,0,loc_vol*sizeof(su3));
  
  //loop over local sites
  for(int i=0;i<loc_vol;i++)
    {
      //Calculate \sum_mu(U_mu(x-mu)-U_mu(x-mu)^dag-U_mu(x)+U^dag_mu(x))
      //and subtract the trace. It is computed just as the TA (traceless anti-herm)
      //this has 8 independent element (A+F+I)=0
      // ( 0,A) ( B,C) (D,E)
      // (-B,C) ( 0,F) (G,H)
      // (-D,E) (-G,H) (0,I)
      for(int mu=0;mu<4;mu++)
	{
	  int b=loclx_neighdw[i][mu];
	  
	  //off-diagonal real parts
	  g[i][0][1][0]+=conf[b][mu][0][1][0]-conf[b][mu][1][0][0]-conf[i][mu][0][1][0]+conf[i][mu][1][0][0]; //B
	  g[i][0][2][0]+=conf[b][mu][0][2][0]-conf[b][mu][2][0][0]-conf[i][mu][0][2][0]+conf[i][mu][2][0][0]; //D
	  g[i][1][2][0]+=conf[b][mu][1][2][0]-conf[b][mu][2][1][0]-conf[i][mu][1][2][0]+conf[i][mu][2][1][0]; //G

	  //off diagonal imag parts
	  g[i][0][1][1]+=conf[b][mu][0][1][1]+conf[b][mu][1][0][1]-conf[i][mu][0][1][1]-conf[i][mu][1][0][1]; //C
	  g[i][0][2][1]+=conf[b][mu][0][2][1]+conf[b][mu][2][0][1]-conf[i][mu][0][2][1]-conf[i][mu][2][0][1]; //E
	  g[i][1][2][1]+=conf[b][mu][1][2][1]+conf[b][mu][2][1][1]-conf[i][mu][1][2][1]-conf[i][mu][2][1][1]; //H

	  //digonal imag parts
	  g[i][0][0][1]+=conf[b][mu][0][0][1]+conf[b][mu][0][0][1]-conf[i][mu][0][0][1]-conf[i][mu][0][0][1]; //A
	  g[i][1][1][1]+=conf[b][mu][1][1][1]+conf[b][mu][1][1][1]-conf[i][mu][1][1][1]-conf[i][mu][1][1][1]; //F
	  g[i][2][2][1]+=conf[b][mu][2][2][1]+conf[b][mu][2][2][1]-conf[i][mu][2][2][1]-conf[i][mu][2][2][1]; //I
	}
      
      //compute the trace
      double T3=(g[i][0][0][1]+g[i][1][1][1]+g[i][2][2][1])/3;
      
      //subtract 1/3 of the trace from each element, so to make traceless the matrix
      g[i][0][0][1]-=T3;
      g[i][1][1][1]-=T3;
      g[i][2][2][1]-=T3;
      
      //fill the other parts
      
      //off-diagonal real parts
      g[i][1][0][0]=-g[i][0][1][0];
      g[i][2][0][0]=-g[i][0][2][0];
      g[i][2][1][0]=-g[i][1][2][0];
      
      //off diagonal imag parts
      g[i][1][0][1]=g[i][0][1][1];
      g[i][2][0][1]=g[i][0][2][1];
      g[i][2][1][1]=g[i][1][2][1];
    }
  
  //fast_fourier_accelerate_fixing(g);

  //multiply by a/2 and add 1
  for(int ivol=0;ivol<loc_vol;ivol++)
    for(int ic1=0;ic1<3;ic1++)
      {
	for(int ic2=0;ic2<3;ic2++)
	  for(int ri=0;ri<2;ri++)
	    g[ivol][ic1][ic2][ri]*=alpha/2;
	g[ivol][ic1][ic1][0]+=1;
      }
      
  //overrelax and unitarize
  overrelax(g,g,3,1.75);
}

//here for future usage
double landau_gauge_fixing_functional(quad_su3 *conf)
{  
  double loc_qual=0;
  double glb_qual;
  
  for(int ivol=0;ivol<loc_vol;ivol++)
    for(int mu=0;mu<4;mu++)
      for(int ic=0;ic<3;ic++)
	loc_qual+=conf[ivol][mu][ic][ic][0];
  
  //global reduction
  MPI_Allreduce(&loc_qual,&glb_qual,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  
  return glb_qual/2/3/4/glb_vol;
}

//compute the quality of the gauge fixing
double landau_gauge_fixing_quality(quad_su3 *conf)
{
  communicate_gauge_borders(conf);

  double loc_omega=0;
  for(int ivol=0;ivol<loc_vol;ivol++)
    {
      su3 delta;
      su3_put_to_zero(delta);
      
      for(int mu=0;mu<4;mu++)
	{
	  int bvol=loclx_neighdw[ivol][mu];

	  //compute delta
	  for(int ic1=0;ic1<3;ic1++)
	    for(int ic2=0;ic2<3;ic2++)
	      {
		delta[ic1][ic2][0]+=
		  +conf[bvol][mu][ic1][ic2][0]
		  -conf[bvol][mu][ic2][ic1][0]
		  -conf[ivol][mu][ic1][ic2][0]
		  +conf[ivol][mu][ic2][ic1][0];
		delta[ic1][ic2][1]+=
		  +conf[bvol][mu][ic1][ic2][1]
		  +conf[bvol][mu][ic2][ic1][1]
		  -conf[ivol][mu][ic1][ic2][1]
		  -conf[ivol][mu][ic2][ic1][1];
	      }

	  //compute trace
	  complex trace={0,0};
	  for(int ic=0;ic<3;ic++)
	    {
	      trace[0]+=delta[ic][ic][0];
	      trace[1]+=delta[ic][ic][1];
	    }
	  trace[0]/=3;
	  trace[1]/=3;
	  
	  //subtract trace
	  for(int ic=0;ic<3;ic++)
	    {
	      delta[ic][ic][0]-=trace[0];
	      delta[ic][ic][1]-=trace[1];
	    }
	}
      
      //compute the trace of the square and summ it to omega
      for(int ic1=0;ic1<3;ic1++)
	for(int ic2=0;ic2<3;ic2++)
	  loc_omega+=
	    +delta[ic1][ic2][0]*delta[ic1][ic2][0]
	    +delta[ic1][ic2][1]*delta[ic1][ic2][1];
    }
  
  //global reduction
  double glb_omega;
  MPI_Allreduce(&loc_omega,&glb_omega,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  
  return glb_omega/glb_vol/3;
}

//perform the landau gauge fixing
void landau_gauge_fixing(quad_su3 *conf_out,quad_su3 *conf_in)
{
  memcpy(conf_out,conf_in,sizeof(quad_su3)*loc_vol);
  
  //fixing transformation
  su3 *fixm=appretto_malloc("fixm",loc_vol+loc_bord,su3);

  //fix iteratively up to reaching required precision
  int iter=0;
  double qual_out=landau_gauge_fixing_quality(conf_out);
  do
    {
      //find the next fixing and compute its quality
      double tin=take_time();
      find_steepest_descent_fixing(fixm,conf_out,alpha);
      gauge_transform_conf(conf_out,fixm,conf_out);
      
      //compute change in quality
      double tint=take_time();
      double qual_in=qual_out;
      qual_out=landau_gauge_fixing_quality(conf_out);      
      double delta_qual=qual_out-qual_in;
      master_printf("Iter %d alpha %lg, quality: %lg (%lg req) delta: %lg; %lg %lg s\n",iter,alpha,qual_out,precision,delta_qual,tint-tin,take_time()-tint);
          
      iter++;
    }
  while(qual_out>=precision);
  
  master_printf("Final quality: %lg\n",qual_out);
  
  appretto_free(fixm);
}

int main(int narg,char **arg)
{
  char conf_in_path[1024];

  //basic mpi initialization
  init_appretto();

  if(narg<2) crash("Use: %s input_file",arg[0]);

  open_input(arg[1]);

  read_str_int("L",&(glb_size[1]));
  read_str_int("T",&(glb_size[0]));
  read_str_str("GaugeConf",conf_in_path,1024);
  read_str_double("Alpha",&alpha);
  read_str_double("Precision",&precision);
  
  close_input();

  //Init the MPI grid 
  init_grid();

  ///////////////////////////////////////////

  quad_su3 *conf=appretto_malloc("Conf",loc_vol+loc_bord,quad_su3);
  quad_su3 *fix_conf=appretto_malloc("Conf2",loc_vol+loc_bord,quad_su3);
  
  read_gauge_conf(conf,conf_in_path);
  communicate_gauge_borders(conf);
  
  landau_gauge_fixing(fix_conf,conf);
  
  master_printf("plaq: %.18g\n",global_plaquette(conf));
  master_printf("plaq: %.18g\n",global_plaquette(fix_conf));

  ///////////////////////////////////////////

  appretto_free(conf);
  appretto_free(fix_conf);
  
  close_appretto();

  return 0;
}
