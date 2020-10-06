#include <math.h>

#include "nissa.hpp"
using namespace std;

#include "../src/types/types.hpp"
#include "../src/types/types_routines.hpp"
#include "../src/operators/Wilson_gluon_Klein_Gordon_operator.hpp"
#include "../src/propagators/Wilson_gluon_propagator.hpp"
#include "../src/propagators/tlSym_gluon_propagator.hpp"

spin1prop *prop;
spin1field *temp1,*temp2;
spin1prop *check_id;

//initialize the program
void init_test(int narg,char **arg)
{
  //Basic mpi initialization
  init_nissa(narg,arg);
  
  //init the grid
  init_grid(8,4);
  
  //allocatepropagators
  prop=nissa_malloc("prop",loc_vol,spin1prop);
  temp1=nissa_malloc("temp1",loc_vol+bord_vol,spin1field);
  temp2=nissa_malloc("temp2",loc_vol,spin1field);
  check_id=nissa_malloc("check_id",loc_vol,spin1prop);
}

//close the program
void close_test()
{
  nissa_free(prop);
  nissa_free(temp1);
  nissa_free(temp2);
  nissa_free(check_id);
  
  close_nissa();
}

void check_id_output()
{
  //compute the point to be printed
  coords ix={1,2,1,2};
  int lx,rx;
  get_loclx_and_rank_of_coord(&lx,&rx,ix);
  
  if(rank==rx)
    {
      printf("\n\nK*G-1 at point: (%d,%d,%d,%d)=%d, rank: %d\n",ix[0],ix[1],ix[2],ix[3],lx,rx);
      spinspin_print(check_id[lx]);
    }
  
  //take the squared norm of check_id
  double loc_d=0;
  NISSA_LOC_VOL_LOOP(imom)
  {
    double point_d=0;
    for(int id=0;id<32;id++)
      {
	double t=((double*)check_id[imom])[id];
	point_d+=t*t;
      }
    loc_d+=point_d;
    if(fabs(point_d)>1.e-10)
      master_printf("%d %d %d %d %d  %lg\n",glb_coord_of_loclx[imom][0],glb_coord_of_loclx[imom][1],glb_coord_of_loclx[imom][2],glb_coord_of_loclx[imom][3],imom,point_d);
  }
  double glb_d=sqrt(glb_reduce_double(loc_d)/glb_vol);
  master_printf("\n\nAverage norm2 of K*G-1: %lg\n\n",glb_d);
}

int main(int narg,char **arg)
{
  init_test(narg,arg);
  
  //anti-periodic boundary condition in one space direction
  double theta[4]={0,0.7,0,0};
  
  //covariant gauge fixing constant
  double alpha=0.3;
  
  //gluon information
  gluon_info gl=create_Wilson_gluon_info(alpha,theta);
  
   /////////////////////////// check Wilson gluon Klein Gordon operator in momentum space //////////////////////
  
  master_printf("\nCheck in x space\n\n");
  
  compute_mom_space_Wilson_gluon_propagator(prop,gl);
  
  for(int nu=0;nu<4;nu++)
    {
      //take index nu of the propagator
      NISSA_LOC_VOL_LOOP(imom)
	for(int mu=0;mu<4;mu++)
	  memcpy(temp1[imom][mu],prop[imom][mu][nu],sizeof(complex));
      
      //apply the KG operator in momentum space
      apply_Wilson_gluon_mom_Klein_Gordon_operator(temp2,temp1,gl);
      
      //put back index nu of the propagator
      NISSA_LOC_VOL_LOOP(imom)
        {
	  for(int mu=0;mu<4;mu++)
	    memcpy(check_id[imom][mu][nu],temp2[imom][mu],sizeof(complex));
	  
	  //subtract id
	  check_id[imom][nu][nu][0]-=1;
	}
    }
  
  check_id_output();
  
  /////////////////////////// check Wilson gluon Klein Gordon operator in momentum space //////////////////////
  
  master_printf("\nCheck in momentum space\n\n");
  
  compute_x_space_Wilson_gluon_propagator_by_fft(prop,gl);
  
  for(int nu=0;nu<4;nu++)
    {
      //take index nu of the propagator
      NISSA_LOC_VOL_LOOP(ivol)
	for(int mu=0;mu<4;mu++)
	  memcpy(temp1[ivol][mu],prop[ivol][mu][nu],sizeof(complex));
      set_borders_invalid(temp1);
      
      //apply the KG operator in momentum space
      apply_Wilson_gluon_x_Klein_Gordon_operator(temp2,temp1,gl);
      
      //put back index nu of the propagator
      NISSA_LOC_VOL_LOOP(ivol)
        {
	  for(int mu=0;mu<4;mu++)
	    memcpy(check_id[ivol][mu][nu],temp2[ivol][mu],sizeof(complex));
	}
      
      //subtract id from x=0
      if(rank==0) check_id[0][nu][nu][0]-=1;
    }
  
  check_id_output();
  
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  
  close_test();
  
  return 0;
}
