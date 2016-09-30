#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "qed_corr.hpp"

namespace nissa
{
  //print
  std::string qed_corr_meas_pars_t::get_str(bool full)
  {
    std::ostringstream os;
    
    os<<"MeasQedCorr\n";
    os<<base_fermionic_meas_t::get_str(full);
    
    return os.str();
  }
}
