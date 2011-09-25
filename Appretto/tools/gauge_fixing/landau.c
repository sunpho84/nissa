#include "appretto.h"

void landau_delta(su3 delta,quad_su3 *conf,int ivol)
{
  su3 U1;

  su3_put_to_zero(delta);
  
  for(int idir=0;idir<4;idir++)
    {
      int bvol=loclx_neighdw[ivol][idir];
      
      //Calculate U_nu(x-nu)+U_nu(x)
      for(int ic1=0;ic1<3;ic1++)
	for(int ic2=0;ic2<3;ic2++)
	  for(int ri=0;ri<2;ri++)
	    U1[ic1][ic2][ri]=conf[bvol][idir][ic1][ic2][ri]+conf[ivol][idir][ic1][ic2][ri];
    }
  
  //Calculate U_nu(x-nu)+U_nu(x)-h.c. and the trace
  complex tr={0,0};
  for(int ic1=0;ic1<3;ic1++)
    {
      for(int ic2=0;ic2<3;ic2++)
	{
	  delta[ic1][ic2][0]=U1[ic1][ic2][0]-U1[ic2][ic1][0];
	  delta[ic1][ic2][1]=U1[ic1][ic2][1]+U1[ic2][ic1][1];
	}	  
      for(int ri=0;ri<2;ri++) tr[ri]+=delta[ic1][ic1][ri];
    }
  
  //subtract the trace
  for(int ic=0;ic<3;ic++) for(int ri=0;ri<2;ri++) delta[ic][ic][ri]-=tr[ri]/3;
}

double landau_gauge_fixing_quality(quad_su3 *conf)
{  
  complex loc_qual={0,0};
  double glb_qual;
  
  for(int ivol=0;ivol<loc_vol;ivol++)
    {
      //su3 delta;
      
      //landau_delta(delta,conf,ivol);

      //calculate the product
      //for(int ic1=0;ic1<3;ic1++)
      //for(int ic2=0;ic2<3;ic2++)
      //complex_summ_the_conj1_prod(loc_qual,delta[ic1][ic2],delta[ic2][ic1]);
      for(int mu=0;mu<4;mu++)
	{
	  int b=loclx_neighdw[ivol][mu];
	  for(int ic=0;ic<3;ic++)
	    loc_qual[0]+=conf[ivol][mu][ic][ic][0]+conf[b][mu][ic][ic][0];
	}
    }
  MPI_Reduce(loc_qual,&glb_qual,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD); 
  
  return -glb_qual/glb_vol/3;
}

void landau_gauge_fixing(quad_su3 *conf)
{
  init_random(10);

  double alpha=0.08;

  set_eo_geometry();

  for(int iter=0;iter<20000;iter++)
    for(int par=0;par<2;par++)
      {
	//communicate_gauge_borders(conf);

	if(par==0)
	  {
	    //calculate original plaquette and quality
	    double qual,old_qual=qual;
	    qual=landau_gauge_fixing_quality(conf);
	    double plaq=global_plaquette(conf);
	    if(rank==0 && iter) printf("Quality: %g %g\n",qual-old_qual,plaq);
	  }
	
	for(int ivol=0;ivol<loc_vol;ivol++)
	  if(loclx_parity[ivol]==par)
	    {
	      su3 g;
	      
	      //if(iter%100!=99)
		{
		  memset(g,0,sizeof(su3));
		  
		  for(int mu=0;mu<4;mu++)
		    {
		      su3 gt,gd;
		      int b=loclx_neighdw[ivol][mu];
		      su3_subt(gt,conf[b][mu],conf[ivol][mu]);
		      unsafe_su3_hermitian(gd,gt);
		      su3_summ(g,gt,gd);
		    }
	      
		  //subtract the trace and multiply by alpha/2
		  complex tr;
		  su3 expon;
		  su3_trace(tr,g);
		  complex_prod_real(tr,tr,1.0/3);
		  su3_subt_complex(g,g,tr);
		  su3_prod_real(expon,g,alpha/2);
	      
		  //calculate the exponential of the matrix
		  su3_copy(g,expon);
		  su3_summ_real(g,g,1);
		}
		/*
	      else
		{
		  for(int i=0;i<3;i++)
		    {
		      for(int j=0;j<3;j++)
			for(int ri=0;ri<2;ri++)
			  g[i][j][ri]=0.2*ran2(ivol)/(iter/1500+1);
		      g[i][i][0]+=1;
		    }
		    }*/
	      su3_unitarize(g,g);
	    
	      //perform the gauge transf
	      for(int mu=0;mu<4;mu++)
		{
		  int b=loclx_neighdw[ivol][mu];
		  su3 te;
		  
		  su3_prod_su3(te,g,conf[ivol][mu]);
		  su3_copy(conf[ivol][mu],te);
		  
		  su3_prod_su3_dag(te,conf[b][mu],g);
		  su3_copy(conf[b][mu],te);
		}
	    }
      }
}

int main(int narg,char **arg)
{
  char filename[1024];

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
  read_str_str("Filename",filename,1024);

  close_input();

  //Init the MPI grid 
  init_grid();

  ///////////////////////////////////////////

  quad_su3 *conf=allocate_quad_su3(loc_vol+loc_bord,"Conf");
  read_gauge_conf(conf,filename);  
  communicate_gauge_borders(conf);
  
  landau_gauge_fixing(conf);
  
  ///////////////////////////////////////////

  close_appretto();

  return 0;
}
