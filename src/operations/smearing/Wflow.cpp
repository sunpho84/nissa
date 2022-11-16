#ifdef HAVE_CONFIG_H
# include "config.hpp"
#endif

#include "base/vectors.hpp"
#include "communicate/edges.hpp"
#include "new_types/su3_op.hpp"
#include "operations/su3_paths/plaquette.hpp"
#include "routines/mpi_routines.hpp"

#include "Wflow.hpp"

namespace nissa
{
  namespace Wflow
  {
    //we add with the new weight the previous one multiplied by the old weight
    void update_arg(quad_su3 *arg,quad_su3 *conf,double dt,const Coords<bool>& dirs,int iter)
    {
      
      communicate_lx_quad_su3_edges(conf);
      
      //Runge-Kutta coefficients
      constexpr std::array<double,3> RK_wn={1.0/4,8.0/9, 3.0/4};
      constexpr std::array<double,3> RK_wo={0,    -17.0/9, -1};
      
      //add the new argument of the exponential to the old one
      NISSA_PARALLEL_LOOP(ivol,0,locVol)
	FOR_ALL_DIRS(mu)
	  if(dirs(mu))
	    {
	      //compute the new contribution
	      su3 staple,temp;
	      su3_put_to_zero(staple);
	      for(int inu=0;inu<NDIM-1;inu++)
		{
		  const Dir nu=perp_dir[mu.nastyConvert()][inu];
		  const LocLxSite& A=ivol,B=loclxNeighup(A,nu),D=loclxNeighdw(A,nu),E=loclxNeighup(D,mu),F=loclxNeighup(A,mu);
		  unsafe_su3_prod_su3(       temp, conf[A.nastyConvert()][nu.nastyConvert()],conf[B.nastyConvert()][mu.nastyConvert()]);
		  su3_summ_the_prod_su3_dag(staple,temp,       conf[F.nastyConvert()][nu.nastyConvert()]);
		  unsafe_su3_dag_prod_su3(temp,    conf[D.nastyConvert()][nu.nastyConvert()],conf[D.nastyConvert()][mu.nastyConvert()]);
		  su3_summ_the_prod_su3(staple,    temp,       conf[E.nastyConvert()][nu.nastyConvert()]);
		}
	      
	      //build Omega
	      su3 omega;
	      unsafe_su3_prod_su3_dag(omega,staple,conf[ivol.nastyConvert()][mu.nastyConvert()]);
	      
	      //compute Q and weight (the minus is there due to original stout)
	      su3 iQ,Q;
	      unsafe_su3_traceless_anti_hermitian_part(iQ,omega);
	      su3_prod_idouble(Q,iQ,-RK_wn[iter]*dt); //putting here the integration time
	      
	      //combine old and new
	      su3_prod_double(arg[ivol.nastyConvert()][mu.nastyConvert()],arg[ivol.nastyConvert()][mu.nastyConvert()],RK_wo[iter]);
	      su3_summassign(arg[ivol.nastyConvert()][mu.nastyConvert()],Q);
	    }
      NISSA_PARALLEL_LOOP_END;
      THREAD_BARRIER();
    }
    
    //update the conf according to exp(i arg) conf_
    void update_conf(quad_su3 *arg,quad_su3 *conf,const Coords<bool>& dirs)
    {
      
      //integrate
      NISSA_PARALLEL_LOOP(ivol,0,locVol)
	FOR_ALL_DIRS(mu)
	  if(dirs(mu))
	    {
	      su3 expiQ;
	      safe_hermitian_exact_i_exponentiate(expiQ,arg[ivol.nastyConvert()][mu.nastyConvert()]);
	      safe_su3_prod_su3(conf[ivol.nastyConvert()][mu.nastyConvert()],expiQ,conf[ivol.nastyConvert()][mu.nastyConvert()]);
	    }
      NISSA_PARALLEL_LOOP_END;
      set_borders_invalid(conf);
    }
  }
  
  //flow for the given time for a dt using 1006.4518 appendix C
  void Wflow_lx_conf(quad_su3* conf,double dt,const Coords<bool>& dirs)
  {
    //storage for staples
    quad_su3 *arg=nissa_malloc("arg",locVol.nastyConvert(),quad_su3);
    vector_reset(arg);
    
    //we write the 4 terms of the Runge Kutta scheme iteratively
    for(int iter=0;iter<3;iter++)
      {
	Wflow::update_arg(arg,conf,dt,dirs,iter);
	Wflow::update_conf(arg,conf,dirs);
      }
    
    nissa_free(arg);
  }
}
