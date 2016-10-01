#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "qed_corr.hpp"

namespace nissa
{
  //measure the quark number and its derivative w.r.t mu
  void measure_quark_rendens(quad_su3 **conf,theory_pars_t &theory_pars,qed_corr_meas_pars_t &meas_pars,int iconf,int conf_created)
  {
    //open the file, allocate point result and source
    FILE *file=open_file(meas_pars.path,conf_created?"w":"a");
    close_file(file);
  }
  
  //print
  std::string qed_corr_meas_pars_t::get_str(bool full)
  {
    std::ostringstream os;
    
    os<<"MeasQedCorr\n";
    os<<base_fermionic_meas_t::get_str(full);
    
    return os.str();
  }
}
