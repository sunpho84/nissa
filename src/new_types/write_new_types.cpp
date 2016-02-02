#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <string.h>

#include "base/debug.hpp"
#include "base/vectors.hpp"
#include "io/input.hpp"
#include "operations/stag/nucleon.hpp"

namespace nissa
{
  int master_fprintf(FILE *stream,const char *format,...);
  
  //print quark content
  std::string quark_content_t::get_str(bool full)
  {
    std::ostringstream os;
    os<<"Quark\t\t=\t\""<<name.c_str()<<"\"\n";
    if(full||deg!=def_deg()) os<<" Degeneracy\t=\t"<<deg<<"\n";
    if(full||mass!=def_mass()) os<<" Mass\t\t=\t"<<mass<<"\n";
    if(full||re_pot!=def_re_pot()) os<<" RePotCh\t=\t"<<re_pot<<"\n";
    if(full||im_pot!=def_im_pot()) os<<" ImPotCh\t=\t"<<re_pot<<"\n";
    if(full||charge!=def_charge()) os<<" ElecCharge\t=\t"<<charge<<"\n";
    
    return os.str();
  }
  
  //print em_field_pars
  std::string em_field_pars_t::get_str(bool full)
  {
    std::ostringstream os;
    if(full||flag||is_nonstandard())
      {
	os<<"BkgrdEMField\n";
	if(full||fabs(E[0])>1e-14) os<<" Ex\t\t=\t"<<E[0]<<"\n";
	if(full||fabs(E[1])>1e-14) os<<" Ey\t\t=\t"<<E[1]<<"\n";
	if(full||fabs(E[2])>1e-14) os<<" Ez\t\t=\t"<<E[2]<<"\n";
	if(full||fabs(B[0])>1e-14) os<<" Bx\t\t=\t"<<B[0]<<"\n";
	if(full||fabs(B[1])>1e-14) os<<" By\t\t=\t"<<B[1]<<"\n";
	if(full||fabs(B[2])>1e-14) os<<" Bz\t\t=\t"<<B[2]<<"\n";
      }
    
    return os.str();
  }
}
