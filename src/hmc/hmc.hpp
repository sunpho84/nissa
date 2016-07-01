#ifndef _HMC_HPP
#define _HMC_HPP

namespace nissa
{
  //Omelyan lambda
  const double omelyan_lambda=0.1931833;
  //Rat approx
  const int nappr_per_quark=3;
  enum rat_approx_hmc_t{RAT_APPR_PF_GEN,RAT_APPR_QUARK_ACTION,RAT_APPR_QUARK_FORCE};
}

#endif
