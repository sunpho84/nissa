#include "nissa.hpp"

using namespace nissa;

int L=24,T=96;

void find_tm_polvect(spin out,double kappa,double mass,momentum_t p,int r)
{
  //dirac operator in the physical basis
  spinspin op;
  spinspin_put_to_diag(op,mass);
  spinspin_dirac_summ_the_prod_idouble(op,base_gamma+5,tau3[r]*(cosh(p[0])+cos(p[1])+cos(p[2])+cos(p[3])-0.5/kappa));
  spinspin_dirac_summ_the_prod_double(op,base_gamma+map_mu[0],-sinh(p[0]));
  for(int mu=1;mu<NDIM;mu++) spinspin_dirac_summ_the_prod_idouble(op,base_gamma+map_mu[mu],sin(p[mu]));
  
  complex d;
  spinspin_det(d,op);
  printf("%lg +I %lg\n",d[0],d[1]);
  spinspin_print(op);
  master_printf("\n");
}

void find_tm_polvect(spin out,tm_quark_info qu,int r)
{
  momentum_t p;
  p[0]=tm_quark_energy(qu,0);
  for(int mu=1;mu<NDIM;mu++) p[mu]=qu.bc[mu]*M_PI/glb_size[mu];
  find_tm_polvect(out,qu.kappa,qu.mass,p,r);
}

double effective_mass(double a,double b,int t,int TH,int par=1)
{
  double targ=b/a;
  
  double yl;
  double yr;
  
  //increment the range up to reaching opposite sign
  bool q;
  double xl=0,xr=1;
  yl=1-targ;
  int niter=0,niter_max=100000;
  do
    {
      if(par==1) yr=cosh((xr)*(TH-(t+1)))/cosh((xr)*(TH-t))-targ;
      else yr=sinh((xr)*(TH-(t+1)))/sinh((xr)*(TH-t))-targ;
      
      q=((yl<0 && yr<0) || (yl>=0 && yr>=0));
      if(q) xr*=2;
      niter++;
    }
  while(q && niter<niter_max);
  if(niter==niter_max) crash("not converging bracketting");
  
  //bisect
  double ym;
  double m;
  niter=0;
  do
    {
      m=(xl+xr)/2;
      if(par==1) ym=cosh(m*(TH-(t+1)))/cosh(m*(TH-t))-targ;
      else       ym=sinh(m*(TH-(t+1)))/sinh(m*(TH-t))-targ;
      if((yl<0 && ym<0) || (yl>0 && ym>0))
        xl=m;
      else
        xr=m;
      niter++;
    }
  while(fabs(ym)>1.e-14 && niter<niter_max);
  if(niter==niter_max) crash("not converging bisecting");
  
  return m;
}

void in_main(int narg,char **arg)
{
  init_grid(T,L);

  tm_quark_info qu;
  qu.r=0;
  qu.mass=0.4;
  qu.kappa=0.11;
  qu.bc[0]=0;
  for(int mu=1;mu<4;mu++) qu.bc[mu]=0.16;
  master_printf("m0: %lg, E: %+016.016lg\n",m0_of_kappa(qu.kappa),tm_quark_energy(qu,0));
  
  spin po;
  find_tm_polvect(po,qu,0);

  //compute lepton prop
  double prop_time=-take_time();
  spinspin *prop=nissa_malloc("prop",loc_vol,spinspin);
  compute_x_space_twisted_propagator_by_fft(prop,qu);
  prop_time+=take_time();
  master_printf("prop time: %lg s\n",prop_time);
  
  double k=qu.bc[1];
  complex c[16][T];
  memset(c,0,sizeof(complex)*16*T);
  NISSA_LOC_VOL_LOOP(ivol)
  {
    //compute the phase
    double ph=0;
    for(int mu=1;mu<4;mu++) ph+=glb_coord_of_loclx[ivol][mu];
    ph*=-k*M_PI/L;
    complex ex={cos(ph),sin(ph)};
    
    for(int ig=0;ig<16;ig++)
      {
	complex tr;
	trace_spinspin_with_dirac(tr,prop[ivol],base_gamma+ig);
	complex_summ_the_prod(c[ig][glb_coord_of_loclx[ivol][0]],tr,ex);
      }
  }
  
  //compute the prop
  double p[T];
  for(int t=0;t<T;t++) p[t]=c[0][t][RE];//+c[4][t][RE];

  //squared prop
  if(1)
    {
      memset(p,0,sizeof(double)*T);
      for(int t=0;t<T;t++)
	{
	  for(int ig=0;ig<16;ig++)
	    p[t]+=squared_complex_norm(c[ig][t]);
	  p[t]=sqrt(p[t]);
	}
    }

  //write the propagator
  FILE *fpr=open_file("/tmp/prop","w");
  for(int t=0;t<T;t++) master_fprintf(fpr,"%+016.016lg\n",p[t]);
  fclose(fpr);
  
  //compute effective mass
  double eff[T/2-1];
  for(int t=0;t<T/2-1;t++) eff[t]=effective_mass(p[t],p[t+1],t,T/2,-1);
  
  FILE *fef=open_file("/tmp/eff","w");
  for(int t=0;t<T/2-1;t++) master_fprintf(fef,"%+016.016lg\n",eff[t]);
  fclose(fef);
  
  nissa_free(prop);
}

int main(int narg,char **arg)
{
  init_nissa_threaded(narg,arg,in_main);
  close_nissa();
  
  return 0;
}
