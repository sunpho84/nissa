#ifndef _WFLOW_HPP
#define _WFLOW_HPP

#include "routines/ios.hpp"

namespace nissa
{
  //structure to Wilson flow
  struct Wflow_pars_t
  {
    double T;
    double dt;
    double def_T(){return 10;}
    double def_dt(){return 0.2;}
    
    int master_fprintf(FILE *fout,bool full=false)
    {
      int nprinted=0;
      
      if(full||is_nonstandard())
	{
	  nprinted+=nissa::master_fprintf(fout,"WFlow\n");
	  if(full||T!=def_T()) nprinted+=nissa::master_fprintf(fout," FlowTime\t=\t%lg\n",T);
	  if(full||dt!=def_dt()) nprinted+=nissa::master_fprintf(fout," InteStep\t=\t%lg\n",dt);
	}
      
      return nprinted;
    }
    
    int is_nonstandard()
    {
      return
	T!=def_T()||
	dt!=def_dt();
    }
    
    Wflow_pars_t() :
      T(def_T()),
      dt(def_dt()) {}
  };
  
  void adaptative_stout_lx_conf(quad_su3 *conf,double *t,double Tmax,double *ext_dt);
  void Wflow_lx_conf(quad_su3 *conf,double dt);
}

#endif
