#include "nissa.hpp"
using namespace std;

#include "../src/types/types.hpp"
#include "../src/types/types_routines.hpp"
#include "../src/diagrams/tadpole.hpp"

int L;

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

int main(int narg,char **arg)
{
  //Basic mpi initialization
  init_nissa(narg,arg);
  
  if(narg<2) CRASH("use: %s input",arg[0]);
  
  open_input(arg[1]);
  
  //lattice size
  int T,L;
  read_str_int("T",&T);
  read_str_int("L",&L);
  
  //quark info
  momentum_t quark_bc;
  read_str_momentum_t("QuarkBc",quark_bc);
  double kappa,mass;
  read_str_double("Kappa",&kappa);
  read_str_double("Mass",&mass);
  
  //gluon info
  momentum_t gluon_bc;
  read_str_momentum_t("GluonBc",gluon_bc);
  double alpha;
  read_str_double("Alpha",&alpha);
  char gluon_type[1024];
  read_str_str("GluonType",gluon_type,1024);
  
  close_input();
  
  //init the grid
  init_grid(T,L);

  //create quark info
  quark_info quark=create_twisted_quark_info(kappa,mass,quark_bc);
  
  //create gluon ino
  gluon_info gluon;
  if(strcasecmp(gluon_type,"tlSym")==0) gluon=create_tlSym_gluon_info(alpha,gluon_bc);
  else
    if(strcasecmp(gluon_type,"Wilson")==0) gluon=create_Wilson_gluon_info(alpha,gluon_bc);
    else CRASH("Unknown gluon type %s",gluon_type);
  
  ////////////////////////////////////
  
  spinspin *tad=nissa_malloc("tad",loc_vol,spinspin);
  
  compute_tadpole_diagram_in_mom_space(tad,quark,gluon);
  //pass_spinspin_from_x_to_mom_space(tad,tad,quark.bc);

  FILE *fout=open_file("tad","w");
  NISSA_LOC_VOL_LOOP(imom)
    {
      //compute momentum
      double a2p2=0;
      for(int mu=0;mu<4;mu++)
	{
	  double ap=M_PI*(2*glb_coord_of_loclx[imom][mu]+quark.bc[mu])/glb_size[mu];
	  a2p2+=ap*ap;
	}
      
      if(a2p2<=4)
	master_fprintf(fout,"%lg\t%lg\n",a2p2,glb_vol*real_part_of_trace_with_id(tad,imom));
    }
  
  if(rank==0) fclose(fout);
  
  ////////////////////////////////////
  
  nissa_free(tad);
  
  close_nissa();
  
  return 0;
}
