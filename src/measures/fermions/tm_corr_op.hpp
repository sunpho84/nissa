#ifndef _TM_CORR_OP_HPP
#define _TM_CORR_OP_HPP

#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "base/vectors.hpp"
#include "geometry/geometry_lx.hpp"
#include "hmc/theory_pars.hpp"
#include "operations/su3_paths/clover_term.hpp"
#include "operations/su3_paths/plaquette.hpp"
#include "inverters/twisted_clover/cg_invert_tmclovD_eoprec.hpp"
#include "inverters/twisted_mass/cg_invert_tmD_eoprec.hpp"
#include "new_types/su3_op.hpp"
#include "threads/threads.hpp"

namespace nissa
{
  /// Structure to hold all info to invert
  struct tm_corr_op
  {
    /// Chekc if we need clover term
    bool need_clov;
    
    /// Configuration
    quad_su3 *conf;
    
    /// Temporary vector
    spincolor* tmp;
    
    /// Clover term
    clover_term_t *Cl;
    
    /// Inverse clover term
    inv_clover_term_t *invCl;
    
    // Residue
    const double residue;
    
    /// Quark content
    theory_pars_t& tp;
    
    /// Current quark
    int cur_flav;
    
    /// Set internal quark for flavor
    void set_for_quark(const int iflav)
    {
      if(cur_flav!=-1 and cur_flav!=iflav)
	{
	  const quark_content_t& cq=tp.quarks[cur_flav];
	  
	  if(cq.cSW)
	    {
	      chromo_operator_remove_cSW(Cl,cq.cSW);
	      master_printf("Remove cSW for flav %d\n",cur_flav);
	    }
	  
	  rem_backfield_without_stagphases_from_conf(conf,tp.backfield[cur_flav]);
	  master_printf("Remove backfield %d\n",cur_flav);
	}
      
      if(cur_flav!=iflav)
	{
	  cur_flav=iflav;
	  
	  const quark_content_t& q=tp.quarks[iflav];
	  add_backfield_without_stagphases_to_conf(conf,tp.backfield[iflav]);
	  master_printf("Adding backfield %d\n",iflav);
	  
	  master_printf("Plaquette: %.16lg\n",global_plaquette_lx_conf(conf));
	  
	  if(q.cSW)
	    {
	      master_printf("Adding cSW for flav %d\n",iflav);
	      chromo_operator_include_cSW(Cl,q.cSW);
	      invert_twisted_clover_term(invCl,q.mass,q.kappa,Cl);
	    }
	}
    }
    
    /// Command to invert
    void inv(spincolor *out,spincolor *in,const int iflav,const int r)
    {
      set_for_quark(iflav);
      
      const quark_content_t& q=tp.quarks[iflav];
      const nissa::dirac_matr& P=
	(tau3[r]==+1)?Pplus:Pminus;
      
      safe_dirac_prod_spincolor(tmp,P,in);
      if(q.cSW) inv_tmclovD_cg_eoprec(out,NULL,conf,q.kappa,Cl,invCl,q.cSW,q.mass*tau3[r],1000000,residue,tmp);
      else inv_tmD_cg_eoprec(out,NULL,conf,q.kappa,q.mass*tau3[r],1000000,residue,tmp);
      safe_dirac_prod_spincolor(out,P,out);
    }
    
    /// Constructor
    tm_corr_op(eo_ptr<quad_su3>& ext_conf,const double& residue,theory_pars_t& tp) :
      residue(residue),tp(tp),cur_flav(-1)
    {
      for(auto& q : tp.quarks)
	if(q.discretiz!=ferm_discretiz::ROOT_TM_CLOV)
	  crash("not defined for non-Wilson quarks");
      
      // Check if clover is actually different from zero
      need_clov=false;
      for(auto& q : tp.quarks)
	need_clov|=(q.cSW!=0);
      
      tmp=nissa_malloc("tmp",locVol+bord_vol,spincolor);
      
      conf=nissa_malloc("conf",locVol+bord_vol+edge_vol,quad_su3);
      paste_eo_parts_into_lx_vector(conf,ext_conf);
      
      if(need_clov)
	{
	  Cl=nissa_malloc("Cl",locVol+bord_vol,clover_term_t);
	  invCl=nissa_malloc("invCl",locVol+bord_vol,inv_clover_term_t);
	  chromo_operator(Cl,conf);
	}
    }
    
    /// Destructor
    ~tm_corr_op()
    {
      nissa_free(conf);
      
      nissa_free(tmp);
      
      if(need_clov)
	{
	  nissa_free(Cl);
	  nissa_free(invCl);
	}
    }
    
    /// Insert a gamma matrix
    static void ins(spincolor *out,const int igamma,spincolor *in);
    
    /// Contract two quarks with a specific gamma
    static void undiluted_meson_contr(complex* contr,
				      spincolor *bw,
				      spincolor *fw,
				      const int& igamma,
				      const int source_coord);
    
    /// Compute the barionic contractions, for a specific projection
    ///
    /// The two gammas refers to the source and sink projection. The two
    /// Wick contractions are stored separately. The three quarks are
    /// intended to be l-d-l' (like-dislike-like), and the two Wick
    /// contractions are the direct and exchange respectively. The
    /// second Wick contraction makes sense only if l' and l are the
    /// same. Propagators are expected to contain a phase in time, and
    /// are automatically shifted to zero time. Each propagator is a
    /// list of 12 spincolor vector, each corresponding to the twelve
    /// components.
    static void compute_baryon_2pts_proj_contr(complex* contr,             ///< Output, with indices {t,iWick}
					       const int& igSo,            ///< Gamma in the source
					       const int& igSi,            ///< Gamma in the sink
					       spincolor** Q1,             ///< l propagator
					       spincolor** Q2,             ///< d propagator
					       spincolor** Q3,             ///< l' propagator
					       const int source_coord,     ///< Source coordinate
					       const double& temporal_bc); ///< Boundary conditon in time
    
    /// Compute the nucleon contractions
    ///
    /// To compute the proton, pass U as l and D as d. For neutron, the
    /// opposite. The two Wick contractions are internally combined. The
    /// propagators are expected to follow the same convention of above.
    static void compute_nucleon_2pts_contr(complex* contr,             ///< Output, with indices {t,iWick}
					   spincolor** Ql,             ///< l propagator
					   spincolor** Qd,             ///< d propagator
					   const int source_coord,     ///< Source coordinate
					   const double& temporal_bc); ///< Boundary conditon in time
  };
}

#endif
