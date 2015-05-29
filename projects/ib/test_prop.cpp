#include "nissa.hpp"

using namespace nissa;

int L=16,T=32;

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

//compute the effective mass
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

void proj_mom_prop(spinspin *zmp,tm_quark_info le,int imom)
{
  GET_THREAD_ID();

  spinspin *prop=nissa_malloc("prop",loc_vol,spinspin);
  compute_mom_space_twisted_propagator(prop,le);
  int dirs[4]={1,0,0,0};
  pass_spinspin_from_mom_to_x_space(prop,prop,dirs,le.bc);
  
  if(IS_MASTER_THREAD)
    {
      for(int t=0;t<loc_size[0];t++)
	{
	  coords c={t+rank_coord[0]*loc_size[0],glb_coord_of_loclx[imom][1],glb_coord_of_loclx[imom][2],glb_coord_of_loclx[imom][3]};
	  int ivol=glblx_of_coord(c);
	  int glb_t=glb_coord_of_loclx[ivol][0];
	  spinspin_copy(zmp[glb_t],prop[ivol]);
	}
      MPI_Allreduce(MPI_IN_PLACE,zmp,sizeof(spinspin)/sizeof(double)*glb_size[0],MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    }
  THREAD_BARRIER();

  nissa_free(prop);
}

void test_proj(tm_quark_info le,int par_apar)
{
  spinspin pr;
  double e=twisted_particle_anti_particle_projector_of_imom(pr,le,0,par_apar);
  master_printf("energy: %lg\n",e);
  spinspin_print(pr);
  master_printf("pr\n");

  spinspin zmp[glb_size[0]];
  proj_mom_prop(zmp,le,0);
  //for(int t=0;t<glb_size[0]*0+1;t++)
  //{
  spinspin_print(zmp[0]);
  master_printf("zmp1\n");
  spinspin_print(zmp[1]);
  master_printf("zmp1\n");
  //spinspin_print(zmp[2]);
  //master_printf("zmp2\n");
  //}
  
  for(int id=0;id<4;id++)
    for(int jd=0;jd<4;jd++)
      for(int ri=0;ri<2;ri++)
	{
	  printf("%d %d %d  ",id,jd,ri);
	  if(fabs(pr[id][jd][ri])>1e-10) printf("%lg\n",zmp[1][id][jd][ri]/pr[id][jd][ri]);
	  else master_printf("\n");
	}
  
  spin wf[2],wf_bar[2];
  for(int s=0;s<2;s++)
    {
      twisted_wavefunction_of_imom(wf[s],le,0,0,s);
      unsafe_spin_prod_dirac(wf_bar[s],wf[s],base_gamma+map_mu[0]);
      spin_conj(wf_bar[s],wf_bar[s]);
      spin_print(wf[s]);
      master_printf("wf[%d]\n",s);
      spin_print(wf_bar[s]);
      master_printf("wf_bar[%d]\n",s);
      spin t;
      spinspin wwa;
      twisted_on_shell_operator_of_imom(wwa,le,0,false,1);
      unsafe_spinspin_prod_spin(t,wwa,wf[s]);
      master_printf("-%d-\n",s);
      spin_print(t);
      master_printf("---\n");
    }
  master_printf("polarisations\n");
  spinspin ww0,ww1,ww,wwa;
  spin_direct_prod(ww0,wf[0],wf_bar[0]);
  spin_direct_prod(ww1,wf[1],wf_bar[1]);
  spinspin_summ(ww,ww0,ww1);
  spinspin_print(ww);
  master_printf("sum_r\n");
  twisted_on_shell_operator_of_imom(wwa,le,0,true,1);
  spinspin_print(wwa);
  master_printf("on-shell-operator\n");
  spinspin diff;
  spinspin_subt(diff,ww,wwa);
  master_printf("norm_diff: %lg\n",spinspin_norm2(diff));
  
  for(int s=0;s<2;s++)
    {
      twisted_wavefunction_of_imom(wf[s],le,0,1,s);
      unsafe_spin_prod_dirac(wf_bar[s],wf[s],base_gamma+map_mu[0]);
      spin_conj(wf_bar[s],wf_bar[s]);
      spin_print(wf[s]);
      master_printf("wf[%d]\n",s);
      spin_print(wf_bar[s]);
      master_printf("wf_bar[%d]\n",s);
      spin t;
      twisted_on_shell_operator_of_imom(wwa,le,0,false,-1);
      unsafe_spinspin_prod_spin(t,wwa,wf[s]);
      master_printf("-a%d-\n",s);
      spin_print(t);
      master_printf("---\n");
    }
  master_printf("apolarisations\n");
  spin_direct_prod(ww0,wf[0],wf_bar[0]);
  spin_direct_prod(ww1,wf[1],wf_bar[1]);
  spinspin_summ(ww,ww0,ww1);
  spinspin_prodassign_double(ww,-1);
  spinspin_print(ww);
  master_printf("asum_r\n");
  twisted_on_shell_operator_of_imom(wwa,le,0,true,-1);
  spinspin_print(wwa);
  master_printf("aon-shell-operator\n");
  spinspin_subt(diff,ww,wwa);
  master_printf("norm_diff: %lg\n",spinspin_norm2(diff));

  double xi=0;
  double A=(2+sqr(le.mass))/2;
  for(int t=0;t<glb_size[0];t++)
    {
      double p=M_PI*(2*t+le.bc[0])/glb_size[0];
      xi+=1/(A-cos(p));
    }
  master_printf("xi: %lg\n",xi);
  
  for(int t=0;t<glb_size[0];t++)
    {
      spinspin re,re2;
      unsafe_spinspin_prod_spinspin(re,pr,zmp[t]);
      complex ou,ou2;
      trace_spinspin(ou,re);
      //master_printf("%d %+016.016lg\n",t,ou[RE]);
      unsafe_spinspin_prod_spinspin(re2,pr,re);
      trace_spinspin(ou2,re2);
      //master_printf("%d %+016.016lg %+016.016lg\n",t,ou[RE]*glb_vol/glb_size[0],ou2[RE]*glb_vol/glb_size[0]/(2*le.mass));
      int tp=glb_size[0]-t;
      int ri=RE;
      
      double att=0;
      int ig=5;
      trace_spinspin_with_dirac(ou,zmp[t],base_gamma+ig);
      switch(ig)
	{
	case 5:
	  ri=IM;
	  att=(1-cosh(e))*(exp(-e*t)+exp(-e*tp))/(1-exp(-e*glb_size[0]))/(2*sinh(e));
	  if(t==0) att=+0.5+(1-cosh(e))*(1+exp(-e*glb_size[0]))/(1-exp(-e*glb_size[0]))/(2*sinh(e));
	  break;
	case 0:
	  ri=RE;
	  att=le.mass*(exp(-e*t)+exp(-e*tp))/(1-exp(-e*glb_size[0]))/(2*sinh(e));
	  break;
	case 4:
	  ri=RE;
	  att=sinh(e)*(exp(-e*t)-exp(-e*tp))/(1-exp(-e*glb_size[0]))/(2*sinh(e));
	  if(t==0) att=0;
	  break;
	default:
	  crash("unknown gamma: %d",ig);
	}
      master_printf("  %d %+016.016lg %+016.016lg\n",t,ou[ri]/4*glb_vol/glb_size[0],att);
      //spinspin_print(re);
      //master_printf("boh\n");
      //safe_dirac_prod_spinspin(re,base_gamma+map_mu[0],re);
      //spinspin_print(re);
      //master_printf("boh2\n\n");
    }
  
  for(int imom=0;imom<loc_vol;imom++)
    {
      //check prop and inverse
      //spinspin a;
      //mom_space_twisted_operator_of_imom(a,le,imom);
      //spinspin b;
      //mom_space_twisted_propagator_of_imom(b,le,imom);
      //unsafe_spinspin_prod_spinspin(c,a,b);
      //spinspin_prodassign_double(c,glb_vol);
      //spinspin c;
      //spinspin_print(c);
      //master_printf("\n");
      
      //spinspin_print(a);
      //master_printf("\n");

      //spinspin_print(b);
      //master_printf("\n");
      
      spinspin pr;
      twisted_particle_projector_of_imom(pr,le,0);
      spinspin_print(pr);
      master_printf("\n\n");
    }
}

void in_main(int narg,char **arg)
{
  init_grid(T,L);
  
  tm_quark_info qu;
  qu.r=0;
  qu.mass=0.4;
  qu.kappa=0.125;
  qu.bc[0]=0;
  for(int mu=1;mu<4;mu++) qu.bc[mu]=0.;
  master_printf("mcrit: %lg, E: %+016.016lg\n",m0_of_kappa(qu.kappa),tm_quark_energy(qu,0));

  test_proj(qu,0);
  
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
