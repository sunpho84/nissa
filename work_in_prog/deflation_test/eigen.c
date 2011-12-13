#include "nissa.h"

int seed;
char conf_path[1024];
double kappa,mu,precision;

double ext_eig_max;
quad_su3 *conf;

double vol_spincolor_norm2(spincolor *s)
{
  double loc_norm2=0,norm2;
  double *ds=(double*)s;
  for(int i=0;i<loc_vol*3*4*2;i++)
    {
      loc_norm2+=(*ds)*(*ds);
      ds++;
    }
  MPI_Allreduce(&loc_norm2,&norm2,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  
  return norm2;
}

double vol_spincolor_norm(spincolor *s)
{return sqrt(vol_spincolor_norm2(s));}

double vol_spincolor_diff_norm2(spincolor *s1,spincolor *s2)
{
  double loc_norm2=0,norm2;
  double *ds1=(double*)s1,*ds2=(double*)s2;
  for(int i=0;i<loc_vol*3*4*2;i++)
    {
      double c=(*ds1)-(*ds2);
      loc_norm2+=c*c;
      ds1++;ds2++;
    }
  MPI_Allreduce(&loc_norm2,&norm2,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  
  return norm2;
}

void vol_spincolor_prod_real(spincolor *out,spincolor *in,double r)
{
  double *din=(double*)in,*dout=(double*)out;
  for(int i=0;i<loc_vol*3*4*2;i++)
    {
      (*dout)=r*(*din);
      dout++;din++;
    }
}

void vol_spincolor_summassign_prod_with_real(spincolor *out,spincolor *in,double r)
{
  double *din=(double*)in,*dout=(double*)out;
  for(int i=0;i<loc_vol*3*4*2;i++)
    {
      (*dout)+=r*(*din);
      dout++;din++;
    }
}

void vol_spincolor_subtassign_the_prod_with_complex(spincolor *a,spincolor *b,complex c)
{
  for(int ivol=0;ivol<loc_vol;ivol++)
    for(int is=0;is<4;is++)
      for(int ic=0;ic<3;ic++)
	complex_subt_the_prod(a[ivol][is][ic],b[ivol][is][ic],c);
}

double vol_spincolor_normalize(spincolor *out,spincolor *in)
{
  double norm=sqrt(vol_spincolor_norm2(in)/loc_vol/12);
  vol_spincolor_prod_real(out,in,1/norm);
  
  return norm;
}

void vol_spincolor_scalar_prod(complex glb,spincolor *a,spincolor *b)
{
  complex loc={0,0};
  for(int ivol=0;ivol<loc_vol;ivol++)
    for(int is=0;is<4;is++)
      for(int ic=0;ic<3;ic++)
	complex_summ_the_conj1_prod(loc,a[ivol][is][ic],b[ivol][is][ic]);
  MPI_Allreduce(loc,glb,2,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
}

void init_program(char *path)
{
  open_input(path);
  
  //init the MPI grid 
  int L,T;
  read_str_int("L",&L);
  read_str_int("T",&T);
  init_grid(T,L);

  //initialize the random generator
  read_str_int("Seed",&seed);
  start_loc_rnd_gen(seed);

  //load the gauge conf
  read_str_str("GaugeConfPath",conf_path,1024);
  conf=nissa_malloc("Conf",loc_vol,quad_su3);
  read_gauge_conf(conf,conf_path);

  //load kappa, mu and preciison
  read_str_double("Kappa",&kappa);
  read_str_double("Mu",&mu);
  read_str_double("Precision",&precision);

  close_input();
}

double find_min_eigen(quad_su3 *conf,double kappa,spincolor *source,double mu,double precision,int find_min)
{
  double eig_max;
  
  spincolor *source=nissa_malloc("source",loc_vol+loc_bord,spincolor);
  generate_undiluted_source(source,RND_Z4,-1);
  vol_spincolor_normalize(source,source);
  
  spincolor *temp=nissa_malloc("temp",loc_vol+loc_bord,spincolor);
  spincolor *result=nissa_malloc("result",loc_vol+loc_bord,spincolor);
  
  int niter=0;
  double diff=precision;
  
  const int each=10;
  do
    {
      apply_Q2(result,source,conf,kappa,mu,temp,NULL,NULL);      
      if(find_min)
	{
	  for(int ivol=0;ivol<loc_vol;ivol++)
	    for(int id=0;id<4;id++)
	      for(int ic=0;ic<3;ic++)
		for(int ri=0;ri<2;ri++)
		  result[ivol][id][ic][ri]=ext_eig_max*source[ivol][id][ic][ri]-result[ivol][id][ic][ri];
	}
      
      if(niter%each) eig_max=vol_spincolor_normalize(source,result);
      else eig_max=vol_spincolor_normalize(temp,result);
      
      if(find_min) eig_max=ext_eig_max-eig_max;
      
      if(niter%each==0)
	{
	  diff=vol_spincolor_diff_norm2(temp,source)/loc_vol/12;
	  spincolor *buf=temp;
	  temp=source;
	  source=buf;
	  
	  master_printf("Iter %d, max eigenvalue: %16.16lg, reached precision: %lg out of %lg\n",niter,eig_max,diff,precision);
	}
      
      niter++;
    }
  while(diff>=precision);
  
  return eig_max;
}

int main(int narg,char **arg)
{
  //basic mpi initialization
  init_nissa();
  
  if(narg<2) crash("Use: %s input_file\n",arg[0]);

  init_program(arg[1]);

  spincolor *eigenvect=nissa_malloc("EigenVect",loc_vol+loc_bord,spincolor);
  spincolor *source=nissa_malloc("Source",loc_vol+loc_bord,spincolor);
  
  ///////////////////////////////////////////
  
  generate_undiluted_source(source,RND_Z4,-1);
  find_minimal_eigen(eigenvect,conf,source,kappa,mu,precision,1);
  
  ///////////////////////////////////////////
  
  nissa_free(eigenvect);
  nissa_free(source);
  
  close_nissa();
  
  return 0;
}
