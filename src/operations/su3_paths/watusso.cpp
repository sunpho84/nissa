#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "base/global_variables.hpp"
#include "base/thread_macros.hpp"
#include "communicate/communicate.hpp"
#include "linalgs/linalgs.hpp"
#include "new_types/su3.hpp"
#include "operations/shift.hpp"
#include "operations/su3_paths/arbitrary.hpp"
#ifdef USE_THREADS
 #include "routines/thread.hpp"
#endif

/*
   ___________
  |           |
  |     _     |      | 
  |____|_|    |     /|\  
  |           |      |        | sizeh
  |___________|      mu       |
   
  |szeh|
 */

namespace nissa
{
  //compute watusso of a fixed size
  //output will be: dist+(2*dmax+1)*rho
  THREADABLE_FUNCTION_8ARG(compute_watusso, complex*,out, double*,small, double*,big, quad_su3*,in, int,size, int,mu, int,nu, int,dmax)
  {
    GET_THREAD_ID();
    
    //compute the big square
    su3 *big=nissa_malloc("big",loc_vol+bord_vol,su3);
    su3 *small=nissa_malloc("small",loc_vol+bord_vol,su3);
    su3 *periscope=nissa_malloc("periscope",loc_vol+bord_vol,su3);
    complex *loc_res=nissa_malloc("loc_res",loc_vol,complex);
    
    int sizeh=size/2;

    //compute the big
    int big_steps[2*5]={
      mu,-sizeh,
      nu,size,
      mu,size,
      nu,-size,
      mu,-(size-sizeh)};
    path_drawing_t b;
    compute_su3_path(&b,big,in,big_steps,5);

    //compute the small+lines
    int small_steps[2*5]={
      nu,sizeh+1,
      mu,1,
      nu,-1,
      mu,-1,
      nu,-sizeh};
    path_drawing_t s;
    compute_su3_path(&s,small,in,small_steps,5);
    
    //compute the periscope
    int irho=0;
    for(int rho=0;rho<NDIM;rho++) //orthogonal dir
      if(rho!=mu&&rho!=nu)
	{
	  for(int orie=-1;orie<=1;orie+=2)
	    {
	      //copy the vector
	      vector_copy(periscope,small);
	      path_drawing_t p=s;
	      
	      for(int d=0;d<=dmax;d++)
		{
		  //trace it
		  NISSA_PARALLEL_LOOP(ivol,0,loc_vol) trace_su3_prod_su3(loc_res[ivol],periscope[ivol],big[ivol]);
		  //wait and collapse
		  THREAD_BARRIER();
		  int iout=dmax+orie*d;
		  complex_vector_glb_collapse(out[iout+(2*dmax+1)*irho],loc_res,loc_vol);
		  
		  //elong if needed
		  if(d!=dmax) elong_su3_path(&p,periscope,in,rho,-orie,true);
		}
	    }
	  //increase the perpendicular dimension
	  irho++;
      }
    
    nissa_free(loc_res);
    nissa_free(periscope);
    nissa_free(big);
    nissa_free(small);
  }
  THREADABLE_FUNCTION_END
}
