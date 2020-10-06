#include <math.h>
#include <stdlib.h>

#include "nissa.hpp"
using namespace std;

#include "../src/propagators/twisted_propagator.hpp"
#include "../src/propagators/tlSym_gluon_propagator.hpp"
#include "../src/propagators/Wilson_gluon_propagator.hpp"
#include "../src/types/types_routines.hpp"
#include "../src/routines/read_and_write.hpp"
#include "../src/routines/shift.hpp"
#include "../src/routines/correlations.hpp"
#include "../src/diagrams/propagator_self_energy.hpp"
#include "../src/diagrams/tadpole.hpp"

double real_part_of_trace_with_igamma(spinspin *q,int imom)
{
  complex tr={0,0};
  
  for(int mu=0;mu<4;mu++)
    {
      int nu=map_mu[mu];
      for(int id=0;id<4;id++)
        {
          complex t;
          int ig=base_gamma[nu].pos[id];
          unsafe_complex_prod(t,base_gamma[nu].entr[id],q[imom][ig][id]);
          complex_summassign(tr,t);
        }
    }
  return tr[1]/16;
}

double real_part_of_trace_with_id(spinspin *q,int imom)
{
  double t=0;
  for(int mu=0;mu<4;mu++)
    t+=q[imom][mu][mu][0];
  
  return t/4;
}

//return the angle and the dist
void compute_dist2_max_coord_angle(int &d2,int &mc,double &angle,coords x)
{
  d2=angle=mc=0;
  for(int mu=0;mu<4;mu++)
    {   
      d2+=x[mu]*x[mu];
      angle+=x[mu];
      if(x[mu]>mc) mc=x[mu];
    }
  angle=acos(angle/sqrt(4*d2))*180/M_PI;
}

int main(int narg,char **arg)
{
  //Basic mpi initialization
  init_nissa(narg,arg);
  
  if(narg<2) crash("use: %s L",arg[0]);
  int L=atoi(arg[1]);
  
  init_grid(2*L,L);
  
  //quark
  double quark_theta[4]={1,0,0,0};
  double kappa=1.0/8;
  double mass=0;
  quark_info qu=create_twisted_quark_info(kappa,mass,quark_theta);
  
  //gluon
  double gluon_theta[4]={0,0,0,0};
  double alpha=0;
  gluon_info gl=create_tlSym_gluon_info(alpha,gluon_theta);
  
  spinspin *self=nissa_malloc("self",loc_vol,spinspin);
  spinspin *tad=nissa_malloc("tad",loc_vol,spinspin);
  
  int comp=1;
  
  if(comp)
    {
      compute_tadpole_diagram_in_x_space(tad,qu,gl);
      pass_spinspin_from_x_to_mom_space(tad,tad,qu.bc);
      compute_self_energy_twisted_diagram_in_x_space(self,qu,gl);
      pass_spinspin_from_x_to_mom_space(self,self,qu.bc);
    }
  
  spinspin *q_prop;
  if(comp)
    {
      q_prop=nissa_malloc("q_prop",loc_vol,spinspin);
      compute_mom_space_twisted_propagator(q_prop,qu);
    }
  
  if(comp)
    NISSA_LOC_VOL_LOOP(imom)
      {
	spinspin temp;
	
	safe_spinspin_prod_spinspin(temp,self[imom],q_prop[imom]);
	safe_spinspin_prod_spinspin(self[imom],q_prop[imom],temp);
	
	safe_spinspin_prod_spinspin(temp,tad[imom],q_prop[imom]);
	safe_spinspin_prod_spinspin(tad[imom],q_prop[imom],temp);
      }
  
  if(comp)
    {
      pass_spinspin_from_mom_to_x_space(self,self,qu.bc);
      pass_spinspin_from_mom_to_x_space(tad,tad,qu.bc);
      pass_spinspin_from_mom_to_x_space(q_prop,q_prop,qu.bc);
    }
  
  corr16 *tad_corr=nissa_malloc("tad_corr",loc_vol,corr16);
  if(comp) compute_all_2pts_qdagq_correlations(tad_corr,q_prop,tad);
  else read_corr16(tad_corr,"/Users/francesco/QCD/LAVORI/X/raw_corrections/24/tad_corr");

  corr16 *self_corr=nissa_malloc("self_corr",loc_vol,corr16);
  if(comp) compute_all_2pts_qdagq_correlations(self_corr,q_prop,self);
  else read_corr16(self_corr,"/Users/francesco/QCD/LAVORI/X/raw_corrections/24/self_corr");

  double norm;
  if(comp) norm=glb_vol2;
  else norm=1;
  
  int igamma=5,REIM=RE;
  int np=0;
  coords x;
  for(x[0]=0;x[0]<=glb_size[0]/2;x[0]++)
    for(x[1]=0;x[1]<=glb_size[1]/2;x[1]++)
      for(x[2]=x[1];x[2]<=glb_size[2]/2;x[2]++)
        for(x[3]=x[2];x[3]<=glb_size[3]/2;x[3]++)
	  if(np<3)
	    {
	      int ivol=glblx_of_coord(x);
	      
	      int d2,max_coord;
	      double angle;
	      
	      compute_dist2_max_coord_angle(d2,max_coord,angle,x);
	      if(angle<=30 && max_coord<=4)
		{
		  master_printf("%d %lg %lg\n",d2,self_corr[ivol][igamma][REIM]*norm,tad_corr[ivol][igamma][REIM]*norm);
		  np++;
		}
	    }
  
  nissa_free(self_corr);
  nissa_free(tad_corr);
  
  if(comp)
    {
      nissa_free(q_prop);
      nissa_free(self);
      nissa_free(tad);  
    }
  
  close_nissa();
  
  return 0;
}
