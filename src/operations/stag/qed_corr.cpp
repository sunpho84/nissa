#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "base/thread_macros.hpp"
#include "base/vectors.hpp"
#include "communicate/borders.hpp"
#include "free_theory/free_theory_types.hpp"
#include "free_theory/tlSym_gauge_propagator.hpp"
#include "geometry/geometry_mix.hpp"
#include "operations/stag/qed_corr.hpp"
#include "routines/mpi_routines.hpp"
#include "routines/thread.hpp"

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
  
  namespace stag
  {
    //hold the map
    struct contr_t
    {
      std::string name;
      color **A;
      color **B;
      contr_t(std::string name,color **A,color **B) : name(name),A(A),B(B) {}
    };
    
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
  }
  
  using namespace stag;
  
  //compute and print
  THREADABLE_FUNCTION_5ARG(measure_qed_corr, quad_su3**,conf, theory_pars_t,theory_pars, qed_corr_meas_pars_t,meas_pars, int,iconf, int,conf_created)
  {
    GET_THREAD_ID();
    
    //open the file, allocate point result and source
    FILE *file=open_file(meas_pars.path,conf_created?"w":"a");
    NEW_FIELD_T(source);
    
    //set photon
    gauge_info photon;
    photon.alpha=FEYNMAN_ALPHA;
    for(int mu=0;mu<NDIM;mu++) photon.bc[mu]=0;
    photon.c1=C1_WILSON;
    photon.zms=UNNO_ALEMANNA;
    
    //compute tadpole
    momentum_t tadpole;
    compute_tadpole(tadpole,photon);
    
    //allocate
    spin1field *photon_field[2]={nissa_malloc("photon_phi_ev",loc_volh+bord_volh,spin1field),nissa_malloc("photon_phi_od",loc_volh+bord_volh,spin1field)};
    NEW_FIELD_T(S);
    NEW_FIELD_T(SS);
    NEW_FIELD_T(AS);
    NEW_FIELD_T(SAS);
    NEW_FIELD_T(ASAS);
    NEW_FIELD_T(SASAS);
    NEW_FIELD_T(TS);
    NEW_FIELD_T(STS);
    
    //write the map
    std::vector<contr_t> contr_map;
    contr_map.push_back(contr_t("00",S,S));
    contr_map.push_back(contr_t("0S",S,SS));
    contr_map.push_back(contr_t("0T",S,STS));
    contr_map.push_back(contr_t("0M",S,SASAS));
    contr_map.push_back(contr_t("LL",SAS,SAS));
    
    //init the contr
    double *glb_contr=nissa_malloc("glb_contr",nthreads*glb_size[0]*contr_map.size(),double);
    double *loc_contr=glb_contr+thread_id*glb_size[0]*contr_map.size();
    
    for(int icopy=0;icopy<meas_pars.ncopies;icopy++)
      {
	get_eo_photon(photon_field,photon);
	fill_source(source,0);
	
	for(int iflav=0;iflav<theory_pars.nflavs();iflav++)
	  {
	    //base
	    MINV(S,iflav,source);
	    //scalar insertion
	    MINV(SS,iflav,S);
	    //photon insertion
	    insert_external_source(AS,conf,&theory_pars,iflav,photon_field,S,-1);
	    MINV(SAS,iflav,AS);
	    //photon insertion
	    insert_external_source(ASAS,conf,&theory_pars,iflav,photon_field,SAS,-1);
	    MINV(SASAS,iflav,ASAS);
	    //tadpole insertion
	    insert_tadpole(TS,conf,&theory_pars,iflav,S,tadpole,-1);
	    MINV(STS,iflav,TS);
	    
	    vector_reset(glb_contr);
	    for(size_t icontr=0;icontr<contr_map.size();icontr++)
	      for(int par=0;par<2;par++)
		NISSA_PARALLEL_LOOP(ieo,0,loc_volh)
		  {
		    int ivol=loclx_of_loceo[par][ieo];
		    int t=glb_coord_of_loclx[ivol][0];
		    for(int ic=0;ic<NCOL;ic++)
		      loc_contr[t+glb_size[0]*icontr]+=real_part_of_complex_scalar_prod(contr_map[icontr].A[par][ieo][ic],contr_map[icontr].B[par][ieo][ic]);
		  }
	  }
	
	//reduce
	glb_threads_reduce_double_vect(loc_contr,glb_size[0]*contr_map.size());
	if(IS_MASTER_THREAD) glb_nodes_reduce_double_vect(glb_contr,glb_size[0]*contr_map.size());
	
	//print
	double norm=1.0/(meas_pars.nhits*glb_spat_vol);
	for(size_t icontr=0;icontr<contr_map.size();icontr++)
	  {
	    master_fprintf(file,"\n%s\n",contr_map[icontr].name.c_str());
	    for(int t=0;t<glb_size[0];t++) master_fprintf(file, "%+16.16lg\n",glb_contr[t+glb_size[0]*icontr]*norm);
	  }
      }
    
    //free
    DELETE_FIELD_T(STS);
    DELETE_FIELD_T(TS);
    DELETE_FIELD_T(SASAS);
    DELETE_FIELD_T(ASAS);
    DELETE_FIELD_T(SAS);
    DELETE_FIELD_T(AS);
    DELETE_FIELD_T(SS);
    DELETE_FIELD_T(S);
    DELETE_FIELD_T(source);
    for(int par=0;par<2;par++) nissa_free(photon_field[par]);
    nissa_free(glb_contr);
    
    close_file(file);
  }
  THREADABLE_FUNCTION_END
}
