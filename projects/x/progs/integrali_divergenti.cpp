#include <math.h>
#include <stdlib.h>

#include "nissa.hpp"
using namespace std;

#include "../src/propagators/twisted_propagator.hpp"
#include "../src/diagrams/tadpole.hpp"
#include "../src/types/types_routines.hpp"
#include "../src/routines/read_and_write.hpp"
#include "../src/routines/correlations.hpp"

complex *prop0,*propg,*qprop;

//initialize the program
void init_calc(int narg,char **arg)
{
  //Basic mpi initialization
  init_nissa(narg,arg);
  
  if(narg<3) crash("use %s L T",arg[0]);
  int T=atoi(arg[1]);
  int L=atoi(arg[2]);
  
  //init the grid
  init_grid(T,L);
  
  //allocate propagators
  prop0=nissa_malloc("prop0",loc_vol,complex);
  propg=nissa_malloc("propg",loc_vol,complex);
  qprop=nissa_malloc("qprop",loc_vol,complex);
}

//close the program
void close_calc()
{
  nissa_free(prop0);
  nissa_free(propg);
  nissa_free(qprop);
  
  close_nissa();
}

int main(int narg,char **arg)
{
  init_calc(narg,arg);
  
  vector_reset(prop0);
  vector_reset(propg);
  
  double m0=0,mg=0;
  NISSA_LOC_VOL_LOOP(ik)
    {
      double Mp=0;
      for(int mu=0;mu<4;mu++)
	{
	  double p=2*M_PI*glb_coord_of_loclx[ik][mu]/glb_size[mu];
	  double ph=p/2;
	  double sinph=sin(ph);
	  Mp+=sinph*sinph;
	}
      
      prop0[ik][RE]=((glblx_of_loclx[ik]!=0||m0!=0)?(1/(4*Mp+m0*m0)):0)/glb_vol;
      propg[ik][RE]=((glblx_of_loclx[ik]!=0||mg!=0)?(1/(4*Mp+mg*mg)):0)/glb_vol;
    }
  
  fft4d(prop0,prop0,1,+1,0);
  fft4d(propg,propg,1,+1,0);

  NISSA_LOC_VOL_LOOP(ik) unsafe_complex_prod(qprop[ik],prop0[ik],propg[ik]);
  fft4d(qprop,qprop,1,-1,1);
  NISSA_LOC_VOL_LOOP(ik) complex_prodassign_double(qprop[ik],glb_vol);
  
  double F0=4.36923,Eg=0.577216;
  FILE *fout=open_file("/tmp/f","w");
  FILE *gout=open_file("/tmp/g","w");
  NISSA_LOC_VOL_LOOP(ik)
  {
    double p2=0,sp2=0;
    int isl=1;
    for(int mu=0;mu<4;mu++)
      {
	double p=2*M_PI*std::min(glb_coord_of_loclx[ik][mu],glb_size[mu]-glb_coord_of_loclx[ik][mu])/glb_size[mu];
	p2+=p*p;
	sp2+=sqr(sin(p));
	if(glb_coord_of_loclx[ik][mu]>=glb_size[mu]/4) isl=0;
      }
    if(isl) 
      master_fprintf(fout,"%lg %lg\n",p2,16*M_PI*M_PI*qprop[ik][RE]);
    if(p2>0) master_fprintf(gout,"%lg %lg\n",p2,2-Eg+F0-log(p2));
    //if(isl) master_fprintf(fout,"%lg %lg\n",sp2,prop0[ik][RE]*glb_vol);
    //if(p2>0) master_fprintf(gout,"%lg %lg\n",p2,1/p2);
  }
  close_file(fout);
  
  close_calc();
  
  return 0;
}
