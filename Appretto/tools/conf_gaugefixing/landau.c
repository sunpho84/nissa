#include "appretto.h"

double alpha;
double precision;

//compute the steepest descent gauge fixing transformation
void find_fixing(su3 *g,quad_su3 *conf,double alpha)
{
  //loop over local sites
  for(int ivol=0;ivol<loc_vol;ivol++)
    {
      su3 delta;
      su3_put_to_zero(delta);
      
      //Calculate \sum_mu(U_mu(x-mu)-U_mu(x-mu)^dag-U_mu(x)+U^dag_mu(x)]
      for(int mu=0;mu<4;mu++)
	{
	  int bvol=loclx_neighdw[ivol][mu];
	  
	  for(int ic1=0;ic1<3;ic1++)
	    for(int ic2=0;ic2<3;ic2++)
	      {
		//real part
		delta[ic1][ic2][0]+=
		  +conf[bvol][mu][ic1][ic2][0] //U_mu(x-mu)
		  -conf[bvol][mu][ic2][ic1][0] //U_mu(x-mu)^dag
		  -conf[ivol][mu][ic1][ic2][0] //U_mu(x)
		  +conf[ivol][mu][ic2][ic1][0];//U_mu(x)^dag
		//imag part
		delta[ic1][ic2][1]+=
		  +conf[bvol][mu][ic1][ic2][1] //U_mu(x-mu)
		  +conf[bvol][mu][ic2][ic1][1] //U_mu(x-mu)^dag
		  -conf[ivol][mu][ic1][ic2][1] //U_mu(x)
		  -conf[ivol][mu][ic2][ic1][1];//U_mu(x)^dag
	      }
	}
      
      //compute the trace/3
      complex trace={0,0};
      for(int ic=0;ic<3;ic++)
	{
	  trace[0]+=delta[ic][ic][0];
	  trace[1]+=delta[ic][ic][1];
	}
      trace[0]/=3;
      trace[1]/=3;
      
      //subtract the trace (divide by 1/3)
      for(int ic=0;ic<3;ic++)
	{
	  delta[ic][ic][0]-=trace[0];
	  delta[ic][ic][1]-=trace[1];
	}
      
      //multiply by a/2 and add 1
      for(int ic1=0;ic1<3;ic1++)
	{
	  for(int ic2=0;ic2<3;ic2++)
	    for(int ri=0;ri<2;ri++)
	      delta[ic1][ic2][ri]*=alpha/2;
	  //add 1
	  delta[ic1][ic1][0]+=1;
	}
      
      //unitarize
      su3_unitarize(g[ivol],delta);
    }
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
	    +delta[ic1][ic2][0]*delta[ic2][ic1][0]
	    +delta[ic1][ic2][1]*delta[ic2][ic1][1];
    }
  
  //global reduction
  double glb_omega;
  MPI_Allreduce(&loc_omega,&glb_omega,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  
  return glb_omega/glb_vol/3;
}

//perform the landau gauge fixing
void landau_gauge_fixing(quad_su3 *conf_out,quad_su3 *conf_in)
{
  //fixing transformation
  su3 *fixm=appretto_malloc("fixm",loc_vol+loc_bord,su3);

  //fix iteratively up to reaching required precision
  int iter=0;
  double qual_out=landau_gauge_fixing_quality(conf_in);
  do
    {
      double qual_in=qual_out;
      
      //find the next fixing and compute its quality
      find_fixing(fixm,conf_in,alpha);
      gauge_transform_conf(conf_out,fixm,conf_in);
      qual_out=landau_gauge_fixing_quality(conf_out);
      
      //compute change in quality
      double delta_qual=qual_out-qual_in;
      master_printf("Iter %d alpha %lg, quality: %lg delta: %lg\n",iter,alpha,qual_in,delta_qual);
      memcpy(conf_in,conf_out,sizeof(quad_su3)*loc_vol);
      
      iter++;
    }
  while(qual_out>=precision);
  
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
  
  ///////////////////////////////////////////

  appretto_free(conf);
  appretto_free(fix_conf);
  
  close_appretto();

  return 0;
}
