#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "minmax_eigenvalues.hpp"
#include "../../operations/remez/remez_algorithm.hpp"
#include "../../new_types/rat_approx.hpp"

#ifdef USE_THREADS
 #include "routines/thread.hpp"
#endif

namespace nissa
{
  //measure minmax_eigenvalues
  void measure_minmax_eigenvalues(quad_su3 **conf,theory_pars_t &theory_pars,minmax_eigenvalues_meas_pars_t &meas_pars,int iconf,int conf_created)
  {
    FILE *file=open_file(meas_pars.path,conf_created?"w":"a");
    
    close_file(file);

		rat_approx_t appr;
		double minimum=0.003, maximum=1.0;
		int num=-1, den=2;
		double tollerance=0.01, minerr=0.002;
		generate_approx(appr,minimum,maximum,num,den,minerr,tollerance);
  }
  
  //print
  std::string minmax_eigenvalues_meas_pars_t::get_str(bool full)
  {
    std::ostringstream os;
    
    os<<"MeasMinMaxEigenval\n";
    os<<base_fermionic_meas_t::get_str(full);
    
    return os.str();
  }
}
