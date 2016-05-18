#include <nissa.hpp>

#include "conf.hpp"
#include "contr.hpp"
#include "pars.hpp"
#include "prop.hpp"

// #include <immintrin.h>

using namespace nissa;

/*
  We follow eq.6.21 of Gattringer, pag 131 and compute the two Wick
  contractions separately, as in
  
  eps_{a,b,c} eps_{a',b',c'} (Cg5)_{al',be'} (Cg5)_{al,be}
  (P+-)_{ga,ga'} S_{be',be}{b',b} (
   S_{al',al}{a',a} S_{ga',ga}{c',c} -
   S_{al',ga}{a',c} S_{ga',al}{c',a})
   
     a',al'---------a,al           a',al'--@   @--a,al
       |             |		    |       \ /    |
     b',be'---------b,be           b',be'---------b,be
       |             |		    |       / \    |
     c',ga'---------c,ga	   c',ga'--@   @--c,ga
     
     insertions are labelled as abc on the source (left) side
   
 */

//init everything
void init_simulation(char *path)
{
  //open input file and read it
  open_input(path);
  read_input_preamble();
  read_use_photon_field();
  read_loc_hadr_curr();
  read_ape_smearing_pars();
  read_gaussian_smearing_pars();
  read_photon_pars();
  read_seed_start_random();
  read_free_theory_flag();
  read_random_gauge_transform();
  read_nsources();
  read_ngauge_conf();
  
  //set how to compute propagators, how to make bars and how to
  //combine the different kind of propagators
  set_inversions();
  set_Cg5();
  set_bar_contr_list();
  set_mes_contr_list();
  mes_gamma_list.push_back(idirac_pair_t(5,5)); //P5P5
  
  //allocate
  conf=nissa_malloc("conf",loc_vol+bord_vol,quad_su3);
  ape_smeared_conf=nissa_malloc("ape_smeared_conf",loc_vol+bord_vol,quad_su3);
  allocate_photon_fields();
  allocate_source();
  allocate_mes_contr();
  allocate_bar_contr();
  allocate_Q_prop();
}

//init a new conf
void start_new_conf()
{
  setup_conf(conf,old_theta,put_theta,conf_path,rnd_gauge_transform,free_theory);
  
  //reset contractions
  vector_reset(mes_contr);
  vector_reset(bar_contr);
}

//handle to discard the source
void skip_conf()
{
  for(int isource=0;isource<nsources;isource++)
    {
      coords coord;
      generate_random_coord(coord);
    }
}

//close deallocating everything
void close()
{
  master_printf("\n");
  master_printf("Inverted %d configurations.\n",nanalyzed_conf);
  master_printf("Total time: %g, of which:\n",tot_prog_time);
  master_printf(" - %02.2f%s to prepare %d photon stochastic propagators (%2.2gs avg)\n",photon_prop_time/tot_prog_time*100,"%",nphoton_prop_tot,photon_prop_time/nphoton_prop_tot);
  master_printf(" - %02.2f%s to prepare %d generalized sources (%2.2gs avg)\n",source_time/tot_prog_time*100,"%",nsource_tot,source_time/nsource_tot);
  master_printf(" - %02.2f%s to perform %d inversions (%2.2gs avg)\n",inv_time/tot_prog_time*100,"%",ninv_tot,inv_time/ninv_tot);
  master_printf("    of which  %02.2f%s for %d cg inversion overhead (%2.2gs avg)\n",cg_inv_over_time/inv_time*100,"%",ninv_tot,cg_inv_over_time/ninv_tot);
  master_printf(" - %02.2f%s to perform %d baryonic contractions (%2.2gs avg)\n",bar_contr_time/tot_prog_time*100,"%",nbar_contr,bar_contr_time/nbar_contr);
  master_printf(" - %02.2f%s to perform %d mesonic contractions (%2.2gs avg)\n",mes_contr_time/tot_prog_time*100,"%",nmes_contr,mes_contr_time/nmes_contr);
  master_printf(" - %02.2f%s to perform %d smearing (%2.2gs avg)\n",smear_oper_time/tot_prog_time*100,"%",nsmear_oper,smear_oper_time/nsmear_oper);
  
  nissa_free(conf);
  nissa_free(ape_smeared_conf);
  free_photon_fields();
  free_source();
  free_Q_prop();
  free_bar_contr();
  free_mes_contr();
}

double M_of_mom(tm_quark_info qu,double sin2_momh)
{
  double M=m0_of_kappa(qu.kappa)+2*sin2_momh;
  return M;
}

double M_of_mom(tm_quark_info qu,int imom)
{
  momentum_t sin_mom;
  double sin2_mom,sin2_momh;
  get_component_of_twisted_propagator_of_imom(sin_mom,sin2_mom,sin2_momh,qu,imom);
  return M_of_mom(qu,sin2_momh);
}

double mom_comp_of_coord(int ip_mu,tm_quark_info qu,int mu)
{return M_PI*(2*ip_mu+qu.bc[mu])/glb_size[mu];}

double mom_comp_of_site(int ip,tm_quark_info qu,int mu)
{return mom_comp_of_coord(glb_coord_of_loclx[ip][mu],qu,mu);}

void sin_mom(momentum_t sin_mom,int imom,tm_quark_info qu)
{for(int mu=0;mu<NDIM;mu++) sin_mom[mu]=sin(mom_comp_of_coord(glb_coord_of_loclx[imom][mu],qu,mu));}

double sin2_mom(int imom,tm_quark_info qu)
{double out=0;for(int mu=0;mu<NDIM;mu++) out+=sqr(sin(mom_comp_of_coord(glb_coord_of_loclx[imom][mu],qu,mu)));return out;}

double den_of_mom(int imom,tm_quark_info qu)
{
  //takes the momenta part
  momentum_t sin_mom;
  double sin2_mom,sin2_momh;
  get_component_of_twisted_propagator_of_imom(sin_mom,sin2_mom,sin2_momh,qu,imom);
  
  //compute M and the denominator
  double M=M_of_mom(qu,sin2_momh);
  double den=sin2_mom+sqr(M)+sqr(qu.mass);
  
  return den;
}

double mom_prod(momentum_t p,momentum_t q)
{
  double pro=0;
  for(int mu=0;mu<NDIM;mu++) pro+=p[mu]*q[mu];
  return pro;
}

void bar_transf(complex *co,tm_quark_info qu)
{
  //for(int p=0;p<glb_size[0];p++)
  //master_printf("%d %.16lg %.16lg\n",p,co[p][0],co[p][1]);
  //master_printf("\n");
  
  tm_quark_info ba=qu;
  ba.bc[0]=QUARK_BOUND_COND*3;
  
  complex c[glb_size[0]];
  for(int x=0;x<glb_size[0];x++)
    {
      complex_put_to_zero(c[x]);
      for(int ip0=0;ip0<glb_size[0];ip0++)
	{
	  double p0=mom_comp_of_coord(ip0,ba,0);
	  complex fact={cos(p0*x),sin(p0*x)};
	  if(x<glb_size[0]/2) complex_summ_the_prod(c[x],co[ip0],fact);
	  else                complex_summ_the_conj1_prod(c[x],co[ip0],fact);
	}
      master_printf("%d %.16lg %.16lg\n",x,c[x][0],c[x][1]);
    }
}

void bar_contr_free(complex *mess,tm_quark_info qu)
{
  double mu2=qu.mass*qu.mass;
  complex co[glb_size[0]];
  memset(co,0,sizeof(complex)*glb_size[0]);
  NISSA_LOC_VOL_LOOP(p)
    {
      double Mp=M_of_mom(qu,p);
      double dp=den_of_mom(p,qu);
      momentum_t sin_p;
      sin_mom(sin_p,p,qu);
      
      NISSA_LOC_VOL_LOOP(q)
        {
	  double Mq=M_of_mom(qu,q);
	  double dq=den_of_mom(q,qu);
	  momentum_t sin_q;
	  sin_mom(sin_q,q,qu);
	  
	  double pq=mom_prod(sin_p,sin_q);
	  double numden=-NCOL*sqr(NDIRAC)*(mu2-Mp*Mq-pq)/(dp*dq*glb_vol*glb_size[0]);
	  
	  for(int r0=0;r0<glb_size[0];r0++)
	    {
	      coords cr;
	      cr[0]=(3*glb_size[0]+r0-glb_coord_of_loclx[p][0]-glb_coord_of_loclx[q][0])%glb_size[0];
	      for(int mu=1;mu<NDIM;mu++) cr[mu]=(3*glb_size[mu]-glb_coord_of_loclx[p][mu]-glb_coord_of_loclx[q][mu])%glb_size[mu];
	      int r=loclx_of_coord(cr);
	      
	      complex_summ_the_prod_double(co[r0],mess[r],numden);
	    }
	}
    }
  
  bar_transf(co,qu);
}

void check_bar()
{
  tm_quark_info qu;
  qu.bc[0]=QUARK_BOUND_COND;
  for(int mu=1;mu<NDIM;mu++) qu.bc[mu]=0;
  qu.kappa=kappa;
  qu.mass=qmass[0];
  qu.r=0;
  
  master_printf(" -----------------baryon direct -------------------- \n");
  
  //compute the propagator traced
  complex *mess=nissa_malloc("mess",glb_vol,complex);
  vector_reset(mess);
  
  NISSA_LOC_VOL_LOOP(r)
    {
      double dr=den_of_mom(r,qu);
      mess[r][RE]=qu.mass/dr/glb_vol;
      mess[r][IM]=-sin(mom_comp_of_site(r,qu,0))/dr/glb_vol;
    }
  
  bar_contr_free(mess,qu);
  
  nissa_free(mess);
}

void check_bar2()
{
  tm_quark_info qu;
  qu.bc[0]=QUARK_BOUND_COND;
  for(int mu=1;mu<NDIM;mu++) qu.bc[mu]=0;
  qu.kappa=kappa;
  qu.mass=qmass[0];
  qu.r=0;
  
  master_printf(" -----------------baryon direct ins on outdiquark -------------------- \n");
  
  //compute the propagator with self-energy insertion, traced
  complex *mess=nissa_malloc("mess",glb_vol,complex);
  vector_reset(mess);
  
  //photon prop
  spin1prop *pho_pro=nissa_malloc("pho_pro",loc_vol,spin1prop);
  NISSA_LOC_VOL_LOOP(p)
    mom_space_tlSym_gauge_propagator_of_imom(pho_pro[p],photon,p);
  
  double mu2=qu.mass*qu.mass;
  NISSA_LOC_VOL_LOOP(p)
    {
      momentum_t sin_p;
      sin_mom(sin_p,p,qu);
      double Mp=M_of_mom(qu,p);
      double pp=mom_prod(sin_p,sin_p),Mpp=Mp*Mp;
      
      NISSA_LOC_VOL_LOOP(q)
        {
	  double Mq=M_of_mom(qu,q);
	  momentum_t sin_q;
	  sin_mom(sin_q,q,qu);
	  
	  //find p-q and compute gluon prop
	  coords cpmq;
	  for(int mu=0;mu<NDIM;mu++) cpmq[mu]=(glb_size[mu]+glb_coord_of_loclx[p][mu]-glb_coord_of_loclx[q][mu])%glb_size[mu];
	  int pmq=glblx_of_coord(cpmq);
	  
	  double pq=mom_prod(sin_p,sin_q);
	  complex num={
	    8*qu.mass*(-Mpp+2*Mp*Mq+mu2-pp+pq),
	    -4*(2*sin_p[0]*(2*(Mp*Mq+mu2)+pq)-(Mpp+mu2+pp)*sin_q[0])};
	  
	  complex_summ_the_prod_double(mess[p],num,pho_pro[pmq][0][0][RE]/den_of_mom(q,qu));
	}
      complex_prodassign_double(mess[p],1/(2*sqr(den_of_mom(p,qu))*glb_vol));
    }
  
  bar_contr_free(mess,qu);
  
  nissa_free(mess);
  nissa_free(pho_pro);
}

void mes_transf(complex *co,tm_quark_info qu)
{
  // for(int p=0;p<glb_size[0];p++)
  //   master_printf("%d %.16lg %.16lg\n",p,co[p][0],co[p][1]);
  // master_printf("\n");
  
  tm_quark_info me=qu;
  me.bc[0]=0;
  
  complex c[glb_size[0]];
  for(int x=0;x<glb_size[0];x++)
    {
      complex_put_to_zero(c[x]);
      for(int ip0=0;ip0<glb_size[0];ip0++)
	{
	  double p0=mom_comp_of_coord(ip0,me,0);
	  complex fact={cos(p0*x),sin(p0*x)};
	  complex_summ_the_prod(c[x],co[ip0],fact);
	}
      master_printf("%d %.16lg %.16lg\n",x,c[x][0],c[x][1]);
    }
}

void check_mes()
{
  tm_quark_info qu;
  qu.bc[0]=QUARK_BOUND_COND;
  for(int mu=1;mu<NDIM;mu++) qu.bc[mu]=0;
  qu.kappa=kappa;
  qu.mass=qmass[0];
  qu.r=0;
  
  master_printf(" ------------------ meson------------------ \n");
  
  double mu2=qu.mass*qu.mass;
  complex co[glb_size[0]];
  memset(co,0,sizeof(complex)*glb_size[0]);
  NISSA_LOC_VOL_LOOP(p)
    {
      momentum_t sin_p;
      sin_mom(sin_p,p,qu);
      double dp=den_of_mom(p,qu);
      double Mp=M_of_mom(qu,p);
      
      for(int q0=0;q0<glb_size[0];q0++)
	{
	  coords cq;
	  cq[0]=q0;
	  for(int mu=1;mu<NDIM;mu++) cq[mu]=glb_coord_of_loclx[p][mu];
	  int q=loclx_of_coord(cq);
	  
	  momentum_t sin_q;
	  sin_mom(sin_q,q,qu);
	  double dq=den_of_mom(q,qu);
	  double Mq=M_of_mom(qu,q);
	  double pq=mom_prod(sin_p,sin_q);
	
	int pmq0=(glb_size[0]+glb_coord_of_loclx[p][0]-glb_coord_of_loclx[q][0])%glb_size[0];
	double contr=4*(mu2+pq+Mp*Mq);
	
	co[pmq0][RE]+=NCOL*contr/(glb_size[0]*glb_vol*dp*dq);
      }
    }
  
  mes_transf(co,qu);
}

void check_mes2()
{
  tm_quark_info qu;
  qu.bc[0]=QUARK_BOUND_COND;
  for(int mu=1;mu<NDIM;mu++) qu.bc[mu]=0;
  qu.kappa=kappa;
  qu.mass=qmass[0];
  qu.r=0;
  
  master_printf(" ------------------ meson with self energy ------------------ \n");
  
  spin1prop *pho_pro=nissa_malloc("pho_pro",loc_vol,spin1prop);
  NISSA_LOC_VOL_LOOP(p)
    mom_space_tlSym_gauge_propagator_of_imom(pho_pro[p],photon,p);
  
  double mu2=qu.mass*qu.mass;
  double mu4=mu2*mu2;
  
  complex co[glb_size[0]];
  memset(co,0,sizeof(complex)*glb_size[0]);
  NISSA_LOC_VOL_LOOP(p)
    {
      momentum_t sin_p;
      sin_mom(sin_p,p,qu);
      double dp=den_of_mom(p,qu);
      double Mp=M_of_mom(qu,p);
      double Mpp=Mp*Mp;
      double pp=mom_prod(sin_p,sin_p);
      
      NISSA_LOC_VOL_LOOP(q)
        {
	  momentum_t sin_q;
	  sin_mom(sin_q,q,qu);
	  double dq=den_of_mom(q,qu);
	  double Mq=M_of_mom(qu,q);
	  double pq=mom_prod(sin_p,sin_q);
	  
	  //find p-q and compute gluon prop
	  coords cpmq;
	  for(int mu=0;mu<NDIM;mu++) cpmq[mu]=(glb_size[mu]+glb_coord_of_loclx[p][mu]-glb_coord_of_loclx[q][mu])%glb_size[mu];
	  int pmq=glblx_of_coord(cpmq);
	  
	  for(int t0=0;t0<glb_size[0];t0++)
	    {
	      coords ct;
	      ct[0]=t0;
	      for(int mu=1;mu<NDIM;mu++) ct[mu]=glb_coord_of_loclx[p][mu];
	      int t=loclx_of_coord(ct);
	      momentum_t sin_t;
	      sin_mom(sin_t,t,qu);
	      double dt=den_of_mom(t,qu);
	      double Mt=M_of_mom(qu,t);
	      double pt=mom_prod(sin_p,sin_t);
	      double qt=mom_prod(sin_q,sin_t);
	      
	      int pmt0=(glb_size[0]+glb_coord_of_loclx[p][0]-t0)%glb_size[0];
	      double num=8*(2*mu4+2*(Mpp*(Mq*Mt-mu2)-Mq*Mt*(mu2+pp)+mu2*(-pp+pq)+(2*mu2+pq)*pt+Mp*(Mt*(2*mu2+pq)+2*Mq*(mu2+pt)))-(Mpp+mu2+pp)*qt);
	      co[pmt0][RE]+=NCOL*num*pho_pro[pmq][0][0][RE]/(dp*dq*dp*dt)/glb_vol/glb_size[0];
	    }
	}
    }
  
  mes_transf(co,qu);
  
  nissa_free(pho_pro);
}

void in_main(int narg,char **arg)
{
  //Basic mpi initialization
  tot_prog_time-=take_time();
  
  //check argument
  if(narg<2) crash("Use: %s input_file",arg[0]);
  
  //init simulation according to input file
  init_simulation(arg[1]);
  
  check_mes();
  master_printf("\n\n");
  check_mes2();
  master_printf("\n\n");
  check_bar();
  master_printf("\n\n");
  check_bar2();
  
  //loop over the configs
  int iconf=0;
  while(read_conf_parameters(iconf,skip_conf,finish_file_present))
    {
      //setup the conf and generate the source
      start_new_conf();
      
      for(int isource=0;isource<nsources;isource++)
	{
	  master_printf("\n=== Source %d/%d ====\n",isource+1,nsources);
	  //shift the conf and create the stochastic photon field
	  random_shift_gauge_conf(conf,old_theta,put_theta);
	  generate_photon_stochastic_propagator();
	  //generate source and smear it
	  generate_original_source();
	  smear_oper_time-=take_time();
	  gaussian_smearing(original_source,original_source,ape_smeared_conf,gaussian_smearing_kappa,gaussian_smearing_niters);
	  smear_oper_time+=take_time();
	  //compute prop and contrelators
	  generate_quark_propagators();
	  compute_bar_contr();
	  compute_mes_contr();
	}
      
      //print out contractions
      print_bar_contr();
      print_mes_contr();
      
      mark_finished();
    }
  
  //close the simulation
  tot_prog_time+=take_time();
  close();
}

int main(int narg,char **arg)
{
  init_nissa_threaded(narg,arg,in_main);
  close_nissa();
  
  return 0;
}
