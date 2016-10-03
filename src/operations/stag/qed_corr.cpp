#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "base/thread_macros.hpp"
#include "base/vectors.hpp"
#include "communicate/borders.hpp"
#include "free_theory/free_theory_types.hpp"
#include "free_theory/tlSym_gauge_propagator.hpp"
#include "geometry/geometry_mix.hpp"
#include "qed_corr.hpp"

namespace nissa
{
  using namespace stag;
  
  //print
  std::string qed_corr_meas_pars_t::get_str(bool full)
  {
    std::ostringstream os;
    
    os<<"MeasQedCorr\n";
    os<<base_fermionic_meas_t::get_str(full);
    
    return os.str();
  }
  
  //return directly a eosplit photon field
  void get_eo_photon(spin1field **out,gauge_info photon)
  {
    //allocate lx version of photon field
    spin1field *photon_eta=nissa_malloc("photon_eta",loc_vol+bord_vol,spin1field);
    spin1field *photon_field=nissa_malloc("photon_field",loc_vol+bord_vol,spin1field);
    spin1field *photon_phi=nissa_malloc("photon_phi",loc_vol+bord_vol,spin1field);
    
    //generate source and stochastich propagator
    generate_stochastic_tlSym_gauge_propagator(photon_phi,photon_eta,photon);
    multiply_by_sqrt_tlSym_gauge_propagator(photon_field,photon_eta,photon);
    split_lx_vector_into_eo_parts(out,photon_field);
    
    nissa_free(photon_phi);
    nissa_free(photon_field);
    nissa_free(photon_eta);
  }
  
  void insert_tadpole_handle(complex out,spin1field **aux,int par,int ieo,int mu,void *pars){out[RE]=((double*)pars)[mu];out[IM]=0;}
  void insert_conserved_current_handle(complex out,spin1field **aux,int par,int ieo,int mu,void *pars){out[RE]=((int*)pars)[mu];out[IM]=0;}
  
  //insert the tadpol
  THREADABLE_FUNCTION_7ARG(insert_tadpole, color**,out, quad_su3**,conf, theory_pars_t*,theory_pars, int,iflav, color**,in, double*,tad, int,t)
  {
    //call with no source insertion, plus between fw and bw, and a global -0.25
    complex fw_factor={-0.25,0},bw_factor={-0.25,0};
    insert_vector_vertex(out,conf,theory_pars,iflav,NULL,in,fw_factor,bw_factor,insert_tadpole_handle,t,tad);
  }
  THREADABLE_FUNCTION_END
  
  //insert the external source, that is one of the two extrema of the stoch prop
  THREADABLE_FUNCTION_7ARG(insert_external_source, color**,out, quad_su3**,conf, theory_pars_t*,theory_pars, int,iflav, spin1field**,curr, color**,in, int,t)
  {
    //call with source insertion, minus between fw and bw, and a global i*0.5
    complex fw_factor={0,+0.5},bw_factor={0,-0.5};
    insert_vector_vertex(out,conf,theory_pars,iflav,curr,in,fw_factor,bw_factor,insert_external_source_handle,t);
  }
  THREADABLE_FUNCTION_END
  
  //compute and print
  void measure_qed_corr(quad_su3 **conf,theory_pars_t &theory_pars,qed_corr_meas_pars_t &meas_pars,int iconf,int conf_created)
  {
    //open the file, allocate point result and source
    FILE *file=open_file(meas_pars.path,conf_created?"w":"a");
    NEW_FIELD_T(source);
    
    //set photon
    gauge_info photon;
    photon.alpha=FEYNMAN_ALPHA;
    for(int mu=0;mu<NDIM;mu++) photon.bc[mu]=0;
    photon.c1=C1_WILSON;
    photon.zms=UNNO_ALEMANNA;
    
    spin1field *photon_field[2]={nissa_malloc("photon_phi_ev",loc_volh+bord_volh,spin1field),nissa_malloc("photon_phi_od",loc_volh+bord_volh,spin1field)};
    NEW_FIELD_T(S);
    NEW_FIELD_T(AS);
    
    for(int icopy=0;icopy<meas_pars.ncopies;icopy++)
      {
	get_eo_photon(photon_field,photon);
	fill_source(source);
	
	for(int iflav=0;iflav<theory_pars.nflavs();iflav++)
	  {
	    MINV(S,iflav,source);
	    insert_external_source(AS,conf,&theory_pars,iflav,photon_field,S,-1);
	  }
      }
    
    nissa_free(photon_field);
    
    close_file(file);
  }
}
