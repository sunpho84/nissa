#include <math.h>
#include <string.h>

#include "nissa.h"

//momentums
momentum_t *mom;
momentum_t *sin_mom;
double *mom2;
double *sin2_mom;

//wilson propagator in momentum and x space
spinspin *mom_space_wi_prop;
spinspin *x_space_wi_prop;

//identical conf
quad_su3 *conf;

//anti-periodic boundary condition in time
double theta[4]={1,0,0,0};

//initialize the program
void init_x()
{
  //Basic mpi initialization
  init_nissa();
  
  //init the grid
  init_grid(8,4);
  
  //allocate momenta
  mom=nissa_malloc("mom",loc_vol,momentum_t); //rothe, 4.7b
  sin_mom=nissa_malloc("sin_mom",loc_vol,momentum_t);
  mom2=nissa_malloc("mom2",loc_vol,double);
  sin2_mom=nissa_malloc("sin2_mom",loc_vol,double);
  
  //allocate wilson propagator in mom and x space
  mom_space_wi_prop=nissa_malloc("mom_space_wi_prop",loc_vol,spinspin);
  x_space_wi_prop=nissa_malloc("mom_space_wi_prop",loc_vol,spinspin);
  
  //identical configuration
  conf=nissa_malloc("id_conf",loc_vol,quad_su3);
  nissa_loc_vol_loop(ivol)
    for(int mu=0;mu<4;mu++)
      su3_put_to_id(conf[ivol][mu]);
  put_boundaries_conditions(conf,theta,0,0);
}

//close the program
void close_x()
{
  nissa_free(mom);
  nissa_free(sin_mom);
  nissa_free(mom2);
  nissa_free(sin2_mom);
  nissa_free(mom_space_wi_prop);
  nissa_free(x_space_wi_prop);
  nissa_free(conf);
  
  close_nissa();
}

//compute the quark propagator in the momentum space
void compute_mom_space_wilson_propagator(spinspin *prop,double kappa)
{
  double M=1/(2*kappa)-4;
  
  //reset the propagator
  memset(prop,0,loc_vol*sizeof(spinspin));
  
  nissa_loc_vol_loop(imom)
    {
      //rothe 4.29b: Mp=M+2*\sum_mu sin^2(p_mu/2)
      double Mp=0;
      for(int mu=0;mu<4;mu++)
	{
	  double ph=mom[imom][mu]/2;
	  double sinph=sin(ph);
	  Mp+=sinph*sinph;
	}
      Mp=2*Mp+M;
      
      //denominator in 4.29a: den=1/(sin2_mom+Mp*Mp)
      double rep_den=1/(sin2_mom[imom]+Mp*Mp);
      
      //compute the ratio using implemented gamma
      for(int ig=0;ig<4;ig++)
	{
	  complex_prod_double(          prop[imom][ig][base_gamma[0].pos[ig]],base_gamma[0].entr[ig],Mp*rep_den);
	  complex_summ_the_prod_idouble(prop[imom][ig][base_gamma[1].pos[ig]],base_gamma[1].entr[ig],-sin_mom[imom][1]*rep_den);
	  complex_summ_the_prod_idouble(prop[imom][ig][base_gamma[2].pos[ig]],base_gamma[2].entr[ig],-sin_mom[imom][2]*rep_den);
	  complex_summ_the_prod_idouble(prop[imom][ig][base_gamma[3].pos[ig]],base_gamma[3].entr[ig],-sin_mom[imom][3]*rep_den);
	  complex_summ_the_prod_idouble(prop[imom][ig][base_gamma[4].pos[ig]],base_gamma[4].entr[ig],-sin_mom[imom][0]*rep_den);
	}
    }
}

//pass from p to x space
void xspace_from_mom_space(complex *out,complex *in,momentum_t theta,int ncps)
{
  //compute the main part of the fft
  fft4d(out,in,ncps,+1,1);

  //compute steps
  momentum_t steps;
  for(int mu=0;mu<4;mu++) steps[mu]=theta[mu]*M_PI/glb_size[mu];

  //add the fractional phase
  nissa_loc_vol_loop(imom)
    {
      //compute local phase
      double arg=0;
      for(int mu=0;mu<4;mu++) arg+=steps[mu]*glb_coord_of_loclx[imom][mu];
      
      //compute the theta
      complex theta={cos(arg),sin(arg)};
      
      //adapt the phase
      for(int ic=0;ic<ncps;ic++)
	safe_complex_prod(out[imom*ncps+ic],out[imom*ncps+ic],theta);
    }  
}

void compute_pion_correlator(complex *corr,spinspin *p)
{
  //stupid propagators
  su3spinspin *stp=nissa_malloc("stp",loc_vol,su3spinspin);
  
  //stupidly promote the propagator to su3
  memset(stp,0,sizeof(su3spinspin)*loc_vol);
  nissa_loc_vol_loop(ivol)
    for(int ic1=0;ic1<1;ic1++)
      for(int ic2=0;ic2<1;ic2++)
	memcpy(stp[ivol][ic1][ic2],p[ivol],sizeof(spinspin));
  
  trace_g_ccss_dag_g_ccss(corr,&(base_gamma[0]),stp,&(base_gamma[0]),stp,1); 
  
  nissa_free(stp);
}

//compute numerically the wilson propagator
void numerically_compute_x_space_wilson_propagator(spinspin *x_space_wi_prop,double kappa)
{
  //source and temp prop
  spincolor *source=nissa_malloc("source",loc_vol,spincolor);
  spincolor *tprop=nissa_malloc("tprop",loc_vol,spincolor);
  
  for(int id_so=0;id_so<4;id_so++)
    {
      //prepare the source
      memset(source,0,sizeof(spincolor)*loc_vol);  
      if(rank==0) source[0][id_so][0][0]=1;
      
      //invert and copy into the spinspin
      inv_tmD_cg_eoprec_eos(tprop,NULL,conf,kappa,0,1000000,1.e-28,source);
      nissa_loc_vol_loop(ivol)
	for(int id_si=0;id_si<4;id_si++)
	  memcpy(x_space_wi_prop[ivol][id_si][id_so],tprop[ivol][id_si][0],sizeof(complex));
    }
  
  //remove the boundary condition
  nissa_loc_vol_loop(ivol)
    {
      //compute the phase
      double arg=0;
      for(int mu=0;mu<4;mu++) arg+=theta[mu]*glb_coord_of_loclx[ivol][mu]/glb_size[mu];
      arg*=M_PI;
      complex phase={cos(arg),sin(arg)};
      
      //multiply for the phase
      for(int id1=0;id1<4;id1++)
	for(int id2=0;id2<4;id2++)
	  safe_complex_prod(x_space_wi_prop[ivol][id1][id2],x_space_wi_prop[ivol][id1][id2],phase);
    }
  
  nissa_free(source);
  nissa_free(tprop);
}

int main(int narg,char **arg)
{
  double kappa=0.177;
  
  init_x();
  
  /////////////////////////////// propagator and pion computed analytically //////////////////////////
  
  define_local_momenta(mom,mom2,sin_mom,sin2_mom,theta);
  
  compute_mom_space_wilson_propagator(mom_space_wi_prop,kappa);
  xspace_from_mom_space((complex*)x_space_wi_prop,(complex*)mom_space_wi_prop,theta,16);
  
  complex pi[glb_size[0]];  
  compute_pion_correlator(pi,x_space_wi_prop);
  int five[1]={5};
  print_contractions_to_file(stdout,1,five,five,pi,0,"",1.0);
  
  print_spinspin(x_space_wi_prop[66]);
  
  /////////////////////////////// propagator and pion computed numerically //////////////////////////
  
  numerically_compute_x_space_wilson_propagator(x_space_wi_prop,kappa);
  compute_pion_correlator(pi,x_space_wi_prop);
  print_contractions_to_file(stdout,1,five,five,pi,0,"",1.0);

  print_spinspin(x_space_wi_prop[66]);
  
  close_x();
  
  return 0;
}
