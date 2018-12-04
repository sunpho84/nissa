#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "minmax_eigenvalues.hpp"
#include "operations/remez/remez_algorithm.hpp"
#include "new_types/rat_approx.hpp"

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
    
    const double minimum=1e-6;
    const double maximum=1.0;
    const int num=-2;
    const int den=3;
    const double tollerance=0.01;
    const double minerr=0.0;
    
    const double t_in=take_time();
    rat_approx_t appr;
    appr.resize(20);
    double res=generate_approx(appr,minimum,maximum,num,den,minerr,tollerance);
    master_printf("Result, res: %lg\n",res);
    appr.master_fprintf(stdout);
    appr.master_fprintf_expr(stdout);
    printf("time required= %.10e secs\n",take_time()-t_in);
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
