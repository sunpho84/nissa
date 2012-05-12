#include <math.h>
#include <string.h>

#include "nissa.h"

//momentums
momentum_t *mom;
momentum_t *sin_mom;
double *mom2;
double *sin2_mom;

//wilson propagator in momentum and x space
spin1prop *mom_space_wi_prop;
spin1prop *x_space_wi_prop;

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
  mom_space_wi_prop=nissa_malloc("mom_space_wi_prop",loc_vol,spin1prop);
  x_space_wi_prop=nissa_malloc("mom_space_wi_prop",loc_vol,spin1prop);
  
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

//compute the wilson action gluon propagator in the momentum space according to P.Weisz
void compute_mom_space_wilson_propagator(spin1prop *prop,double alpha,momentum_t bc)
{
  //check absence of zero modes
  int zmp=1;
  for(int mu=0;mu<4;mu++) zmp&=(bc[mu]==0);
  if(zmp) crash("zero mode present, prop not defined");
  
  //reset the propagator
  memset(prop,0,loc_vol*sizeof(spin1prop));
  
  //delta_mu_nu and A_wilson=1-delta_mu_nu (eq A.5)
  int kron_delta[4][4]={{1,0,0,0},{0,1,0,0},{0,0,1,0},{0,0,0,1}};
  int A[4][4]=         {{0,1,1,1},{1,0,1,1},{1,1,0,1},{1,1,1,0}};
  
  //eq A.4
  nissa_loc_vol_loop(imom)
    {
      //momentum
      momentum_t k,kt;
      double kt2=0;
      for(int mu=0;mu<4;mu++)
	{
	  k[mu]=M_PI*(2*glb_coord_of_loclx[imom][mu]+bc[mu])/glb_size[mu]; //lattice momentum
	  kt[mu]=2*sin(k[mu]/2); //eq 3.7 (factor 2 and 1/2 coming from symmetrized derivative)
	  kt2+=kt[mu]*kt[mu];
	}
      
      for(int nu=0;nu<4;nu++)
	for(int ta=0;ta<4;ta++)
	  {
	    //first term, alpha * kt_nu * kt_ta
	    double p=alpha*kt[nu]*kt[ta];
	    
	    //second term, \sum_si
	    for(int si=0;si<4;si++)
		p+=(kt[si]*kron_delta[ta][nu]-kt[ta]*kron_delta[si][nu])*A[ta][si]*kt[si];
	  
	    //divide everything by kt^2
	    prop[imom][nu][ta][RE]=p/kt2;
	    //imag part is zero
	    prop[imom][nu][ta][IM]=0;
	  }
    }
}

/* klein gordon operator for wilson action, as in P.Weisz eq A.3:
  out(i;mu)=K(i,j;mu,nu)*in(j;nu)
 
 with K being the Klein Gordon operator:
  K(i,j;mu,nu)=Sum_l{D(i,l;mu)*D(l,j;nu)/alpha+Sum_rh[D(i,l;rh)*kron(mu,nu)-D(i,l;nu)*kron(rh,nu)]q(mu,rh)*D(l,j;rh)}

 where we have defined symmetric derivative operator D(i,j;mu)=[kron(j,i+mu)-kron(j,i-mu)]/2
 and we defined: q_mu_nu=1-kron(mu,nu) (this will be different for tree level Symanzik and Iwasaki)

 since this computation involves 2 neighbours it needs to be performed in 2 steps, in particular we will define
  d1(i,mu,nu)=D(i,j;mu)*in(j,nu)

 from which we can compute
  out(i;mu)=D(i,j;mu)*[Sum_nu d1(j;nu,nu)]/alpha+Sum_rh_nu[D(i,j;rh)*kron(mu,nu)-D(i,j;nu)*kron(rh,nu)]*q(mu,rh)*d1(j;rh,nu)
           =D(i,j;mu)*[Sum_nu d1(j;nu,nu)]/alpha+Sum_nu{D(i,j;nu)*q(mu,nu)*[d1(j;nu,mu)-d1(j;nu,nu)]}
*/
void klein_gordon_wilson_action_operator(spin1field *out,spin1field *in,double alpha,momentum_t bc)
{
  spin1prop *d1=nissa_malloc("d1",loc_vol+bord_vol,spin1prop);
  
  //compute the derivative d1
  nissa_loc_vol_loop(ivol)
    for(int mu=0;mu<4;mu++)
      {
	int up_mu=loclx_neighup[ivol][mu];
	int dw_mu=loclx_neighdw[ivol][mu];
	for(int nu=0;nu<4;nu++)
	  {
	    complex_subt(d1[ivol][mu][nu],in[up_mu][nu],in[dw_mu][nu]);
	    complex_prodassign_double(d1[ivol][mu][nu],0.5);
	  }
      }
  
  //border to be communicated
  
  //second step of the computation
  nissa_loc_vol_loop(ivol)
    for(int mu=0;mu<4;mu++)
      {
	int up_mu=loclx_neighup[ivol][mu];
	int dw_mu=loclx_neighdw[ivol][mu];
	
	//write gf term k_mu * k_nu /alpha
	out[ivol][mu][RE]=out[ivol][mu][IM]=0;
	//we must take diagonal part of the derivative, D
	for(int nu=0;nu<4;nu++)
	  {
	    complex_summassign(out[ivol][mu],d1[up_mu][nu]);
	    complex_subtassign(out[ivol][mu],d1[dw_mu][nu]);
	  }
	complex_prodassign_double(out[ivol][mu],0.5/alpha);
	
	//now if q[mu][nu]!=0, i.e. mu!=nu, add the other piece
	for(int nu=0;nu<4;nu++)
	  {
	    complex_summassign(out[ivol][mu],d1[up_mu][nu]);
	    complex_subtassign(out[ivol][mu],d1[dw_mu][nu]);
	  }
	complex_prodassign_double(out[ivol][mu],0.5/alpha);	
      }
  
  nissa_free(d1);
}

//compute the x-space wilson action gluon propagator by inverting the klein-gordon operator
void compute_x_space_wilson_propagator(spin1prop *prop,double alpha,momentum_t bc)
{
  
}

/*
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
*/

int main(int narg,char **arg)
{
  double alpha=0.4;
  
  init_x();
  
  /////////////////////////////// propagator and pion computed analytically //////////////////////////
  
  define_local_momenta(mom,mom2,sin_mom,sin2_mom,theta);
  
  compute_mom_space_wilson_propagator(mom_space_wi_prop,alpha,theta);
  //xspace_from_mom_space((complex*)x_space_wi_prop,(complex*)mom_space_wi_prop,theta,16);
  
  /*
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
  */

  close_x();
  
  return 0;
}
