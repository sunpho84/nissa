#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "base/random.hpp"
#include "base/vectors.hpp"
#include "geometry/geometry_mix.hpp"
#include "hmc/quark_pars.hpp"
#include "inverters/twisted_clover/cg_invert_tmclovD_eoprec.hpp"
#include "inverters/twisted_mass/cg_invert_tmD_eoprec.hpp"
#include "operations/su3_paths/clover_term.hpp"
#include "routines/mpi_routines.hpp"

#include "tm_tuning.hpp"

namespace nissa
{
  namespace
  {
    const int ncorr_kind=4;
  }
  
  //compute correlation functions for twisted clover, needed to fix tuning
  THREADABLE_FUNCTION_6ARG(tm_tuning, complex*,corr, quad_su3*,conf, clover_term_t*,Cl, inv_clover_term_t*,invCl, quark_content_t*,q, tm_tuning_meas_pars_t*,meas_pars)
  {
    GET_THREAD_ID();
    
    spincolor *eta=nissa_malloc("eta",loc_vol+bord_vol,spincolor);
    spincolor *phi=nissa_malloc("phi",loc_vol+bord_vol,spincolor);
    spincolor *phi_ins_S=nissa_malloc("phi_ins_S",loc_vol+bord_vol,spincolor);
    spincolor *phi_ins_P=nissa_malloc("phi_ins_P",loc_vol+bord_vol,spincolor);
    
    //Place where to store local result
    complex *loc_corr=new complex[ncorr_kind*glb_size[0]];
    memset(loc_corr,0,sizeof(complex)*ncorr_kind*glb_size[0]);
    
    //Source time
    int tso;
    
    // Command to invert
    auto inv=[&](spincolor *out,spincolor *in)
	     {
	       if(q->cSW) inv_tmclovD_cg_eoprec(out,NULL,conf,q->kappa,Cl,invCl,q->cSW,q->mass,1000000,meas_pars->residue,in);
	       else inv_tmD_cg_eoprec(out,NULL,conf,q->kappa,q->mass,1000000,meas_pars->residue,in);
	     };
    
    // Command to insert
    auto ins=[&](spincolor *out,int igamma,spincolor *in)
	     {
		 NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
		   {
		     unsafe_dirac_prod_spincolor(out[ivol],base_gamma+igamma,in[ivol]);
		   }
		 NISSA_PARALLEL_LOOP_END;
		 
		 set_borders_invalid(out);
	     };
    
    // Command to contract
    auto contr=[&](spincolor *bw,spincolor *fw,int igamma,int icorr)
	       {
		 dirac_matr g=base_gamma[igamma]*base_gamma[5];
		 
		 NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
		      {
			int t=(glb_coord_of_loclx[ivol][0]-tso+glb_size[0])%glb_size[0];
			for(int id=0;id<NDIRAC;id++)
			  {
			    complex s{0.0,0.0};
			    for(int ic=0;ic<NCOL;ic++)
			      complex_summ_the_conj1_prod(s,bw[ivol][g.pos[id]][ic],fw[ivol][id][ic]);
			    complex_summ_the_prod(loc_corr[icorr+ncorr_kind*t],s,g.entr[id]);
			  }
		      }
		 NISSA_PARALLEL_LOOP_END;
		 
		 THREAD_BARRIER();
	       };
    
    int nhits=meas_pars->nhits;
    for(int hit=0;hit<nhits;hit++)
      {
	//Random source
	coords coord;
	generate_random_coord(coord);
	
	//Source time
	tso=coord[0];
	generate_undiluted_source(eta,meas_pars->rnd_type,tso);
	
	inv(phi,eta);
	ins(phi_ins_P,5,phi);
	inv(phi_ins_P,phi_ins_P);
	inv(phi_ins_S,phi);
	
	contr(phi,phi,5,0);
	contr(phi,phi,4,1);
	contr(phi,phi_ins_S,4,2);
	contr(phi,phi_ins_P,4,3);
      }
    
    glb_threads_reduce_double_vect((double*)loc_corr,2*ncorr_kind*glb_size[0]);
    if(IS_MASTER_THREAD) glb_nodes_reduce_complex_vect(corr,loc_corr,ncorr_kind*glb_size[0]);
    
    delete [] loc_corr;
    
    nissa_free(eta);
    nissa_free(phi);
    nissa_free(phi_ins_S);
    nissa_free(phi_ins_P);
  }
  THREADABLE_FUNCTION_END
  
  //compute and print
  void measure_tm_tuning(quad_su3 **ext_conf,theory_pars_t &tp,tm_tuning_meas_pars_t &meas_pars,int iconf,int conf_created)
  {
    FILE *file=open_file(meas_pars.path,conf_created?"w":"a");
    
    quad_su3 *conf=nissa_malloc("conf",loc_vol+bord_vol+edge_vol,quad_su3);
    paste_eo_parts_into_lx_vector(conf,ext_conf);
    
    // Check if clover is actually different from zero
    bool need_clov=false;
    for(auto& q : tp.quarks)
      need_clov|=(q.cSW!=0);
    
    clover_term_t *Cl=nullptr;
    inv_clover_term_t *invCl=nullptr;
    if(need_clov)
      {
	Cl=nissa_malloc("Cl",loc_vol+bord_vol,clover_term_t);
	invCl=nissa_malloc("invCl",loc_vol+bord_vol,inv_clover_term_t);
	chromo_operator(Cl,conf);
      }
    
    complex *corr=nissa_malloc("corr",glb_size[0]*ncorr_kind,complex);
    
    int ncopies=meas_pars.ncopies;
    for(int icopy=0;icopy<ncopies;icopy++)
      for(int iflav=0;iflav<tp.nflavs();iflav++)
	{
	  master_fprintf(file," # conf %d ; flv = %d , m = %lg\n",
			 iconf,iflav,tp.quarks[iflav].mass);
	  
	  quark_content_t& q=tp.quarks[iflav];
	  
	  if(q.discretiz!=ferm_discretiz::ROOT_TM_CLOV) crash("not defined for non-Wilson quarks");
	  
	  add_backfield_without_stagphases_to_conf(conf,tp.backfield[iflav]);
	  if(q.cSW)
	    {
	      chromo_operator_include_cSW(Cl,q.cSW);
	      invert_twisted_clover_term(invCl,q.mass,q.kappa,Cl);
	    }
	  
	  verbosity_lv2_master_printf("Evaluating tm tuning for flavor %d/%d, ncopies %d/%d\n",
					  iflav+1,tp.nflavs(),icopy+1,ncopies);
	  
	  tm_tuning(corr,conf,Cl,invCl,&q,&meas_pars);
	  
	  if(q.cSW) chromo_operator_remove_cSW(Cl,q.cSW);
	  rem_backfield_without_stagphases_from_conf(conf,tp.backfield[iflav]);
	  
	  //output
	  for(int t=0;t<glb_size[0];t++)
	    {
	      master_fprintf(file,"%d\t",t);
	      for(int ic=0;ic<ncorr_kind;ic++)
		{
		  complex c;
		  complex_prod_double(c,corr[ic+ncorr_kind*t],1.0/(meas_pars.nhits*glb_spat_vol));
		  master_fprintf(file,"\t%+.016lg , %+.016lg",c[RE],c[IM]);
		}
	      master_fprintf(file,"\n");
	    }
	  
	  master_fprintf(file,"\n");
	}
    
    nissa_free(corr);
    
    nissa_free(conf);
    if(need_clov)
      {
	nissa_free(Cl);
	nissa_free(invCl);
      }
    close_file(file);
  }
  
  //nucleon correlators
  std::string tm_tuning_meas_pars_t::get_str(bool full)
  {
    std::ostringstream os;
    
    os<<"MeasTmTuning\n";
    os<<base_fermionic_meas_t::get_str(full);
    
    return os.str();
  }
}
