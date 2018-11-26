#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <algorithm>

#include "base/thread_macros.hpp"
#include "base/vectors.hpp"
#include "communicate/borders.hpp"
#include "communicate/edges.hpp"
#include "geometry/geometry_mix.hpp"
#include "io/input.hpp"
#include "linalgs/linalgs.hpp"
#include "new_types/complex.hpp"
#include "new_types/float_128.hpp"
#include "new_types/spin.hpp"
#include "new_types/su3.hpp"
#include "operations/fft.hpp"
#include "operations/gaugeconf.hpp"
#include "operations/remap_vector.hpp"
#include "operations/su3_paths/gauge_sweeper.hpp"
#include "operations/smearing/stout.hpp"
#include "operations/smearing/Wflow.hpp"
#include "operations/su3_paths/plaquette.hpp"
#include "routines/ios.hpp"
#include "routines/mpi_routines.hpp"
#include "dirac_operators/tmQ/dirac_operator_tmQ.hpp"
#include "dirac_operators/tmQ2/dirac_operator_tmQ2.hpp"
#include "dirac_operators/stD/dirac_operator_stD.hpp"
#include "eigenvalues/eigenvalues_all.hpp"
#include "eigenvalues/eigenvalues.hpp"
#include "eigenvalues/eigenvalues_autarchic.hpp"
#include "inverters/twisted_mass/cg_invert_tmQ2.hpp"

//FIXME clean includes
#include "base/random.hpp"
#include "base/vectors.hpp"
#include "operations/fermions/stag.hpp"
#include "communicate/borders.hpp"
#include "geometry/geometry_eo.hpp"
#include "hmc/backfield.hpp"
#include "inverters/staggered/cg_invert_stD.hpp"
#include "linalgs/linalgs.hpp"
#include "new_types/su3.hpp"
#include "routines/mpi_routines.hpp"

#ifdef USE_THREADS
 #include "routines/thread.hpp"
#endif

#include "spectral_projectors.hpp"

namespace nissa{
  // This measure will compute the first 'n' eigenvalues (parameter) 
  // and eigenvectors of DD^+ in staggered formulation, in order to build an estimate
  // of the topological charge as Q = \sum_i \bar(u)_i gamma5 u_i,
  // where {u_i} are the first n eigenvectors.

  THREADABLE_FUNCTION_6ARG(fill_eigenpart, complex**,eigvec, quad_su3**,conf,complex*, charge_cut, complex*, DD_eig_val,int,neigs, double,eig_precision)
  {
//    GET_THREAD_ID();
    master_fprintf(stdout,"neigs = %d, eig_precision = %.2e\n", neigs, eig_precision);


//    //wrap the application of DD^+ into an object that can be passed to the eigenfinder
    color *out_tmp_eo[2]={nissa_malloc("out_tmp_EVN",loc_volh+bord_volh,color),nissa_malloc("out_tmp_ODD",loc_volh+bord_volh,color)};
    color *in_tmp_eo[2]={nissa_malloc("in_tmp_EVN",loc_volh+bord_volh,color),nissa_malloc("in_tmp_ODD",loc_volh+bord_volh,color)};
    color *fill_tmp_eo[2]={nissa_malloc("fill_tmp_EVN",loc_volh+bord_volh,color),nissa_malloc("fill_tmp_ODD",loc_volh+bord_volh,color)};
	  vector_reset(out_tmp_eo[EVN]);
    vector_reset(out_tmp_eo[ODD]);
	  vector_reset(in_tmp_eo[EVN]);
	  vector_reset(in_tmp_eo[ODD]);
	  vector_reset(fill_tmp_eo[EVN]);
	  vector_reset(fill_tmp_eo[ODD]);

		quad_u1 *u1b[2];
    u1b[0]=nissa_malloc("u1b",loc_volh+bord_volh,quad_u1);
    u1b[1]=nissa_malloc("u1b",loc_volh+bord_volh,quad_u1);
    init_backfield_to_id(u1b);
    const auto imp_mat=[conf,&in_tmp_eo,out_tmp_eo](complex *out_lx,complex *in_lx){
      split_lx_vector_into_eo_parts(in_tmp_eo,(color*)in_lx);

      evn_apply_stD(out_tmp_eo[EVN],conf,0.01,(color**)in_tmp_eo);
      odd_apply_stD(out_tmp_eo[ODD],conf,0.01,(color**)in_tmp_eo);
      
      evn_apply_stD_dag(in_tmp_eo[EVN],conf,0.01,(color**)out_tmp_eo);
      odd_apply_stD_dag(in_tmp_eo[ODD],conf,0.01,(color**)out_tmp_eo);
      
      paste_eo_parts_into_lx_vector((color*)out_lx,in_tmp_eo);
				
        set_borders_invalid(out_lx);
      };
  


    //parameters of the eigensolver
    const bool min_max=0;
    const int mat_size=loc_vol*NCOL;
    const int mat_size_to_allocate=(loc_vol+bord_vol)*NCOL;
    const int niter_max=100000;
    master_printf("mat_size=%d, mat_size_to_allocate=%d\n",mat_size,mat_size_to_allocate); 
    //wrap the generation of the test vector into an object that can be passed to the eigenfinder
    const auto filler=[&fill_tmp_eo](complex *out_lx)
      {
        generate_fully_undiluted_eo_source(fill_tmp_eo,RND_GAUSS,-1,0);
        paste_eo_parts_into_lx_vector((color*)out_lx,fill_tmp_eo);
      };
    
    //launch the eigenfinder
    double eig_time=-take_time();

 

    eigenvalues_of_hermatr_find(eigvec,DD_eig_val,neigs,min_max,mat_size,mat_size_to_allocate,imp_mat,eig_precision,niter_max,filler);

//    double norm_cut[neigs];
    complex *point_vec = nissa_malloc("point_vec",(loc_vol+bord_vol)*NCOL,complex);
    master_printf("\n\nEigenvalues of D^+D:\n");
    for(int ieig=0; ieig<neigs; ++ieig){
      master_printf("%d (%.16lg,%.16lg)\n",ieig,DD_eig_val[ieig][RE],DD_eig_val[ieig][IM]);
      
      // compute partial sum terms of tr(g5)     
      // save 'eigvec[ieigs]' in staggered format (i.e., 'fill_tmp_eo'),
      // then multiply it with gamma5 and save the result in 
      // 'out_tmp_eo'. The term contributing to the partial sum of tr(g5)
      // will be stored in 'charge_cut[ieig]' as the hermitian product 
      // between 'fill_tmp_eo' and 'out_tmp_eo'. 

      split_lx_vector_into_eo_parts((color**)fill_tmp_eo,(color*)eigvec[ieig]);

      //multiply by gamma5
      apply_stag_op(out_tmp_eo,conf,u1b,stag::GAMMA_5,stag::IDENTITY,fill_tmp_eo);

      //take hermitian product
      stag::summ_the_trace((double*)&charge_cut[ieig],(complex*)point_vec,out_tmp_eo,fill_tmp_eo);
//      stag::summ_the_trace((double*)&norm_cut[ieig],point_vec,fill_tmp_eo,fill_tmp_eo);
//      master_printf("charge_cut[%d] = %.10f, norm_cut[%d] = %.10f\n",ieig,charge_cut[ieig],ieig,norm_cut[ieig]);

    }
    master_printf("\n\n\n");

    eig_time+=take_time();
    master_printf("Eigenvalues time: %lg\n",eig_time);
    
    nissa_free(out_tmp_eo[EVN]);
    nissa_free(out_tmp_eo[ODD]);
    nissa_free(in_tmp_eo[EVN]);
    nissa_free(in_tmp_eo[ODD]);
    nissa_free(fill_tmp_eo[EVN]);
    nissa_free(fill_tmp_eo[ODD]);
    nissa_free(u1b[0]);
    nissa_free(u1b[1]);
    nissa_free(point_vec);
  }
  THREADABLE_FUNCTION_END

  
  //measure the topological charge
  //TODO: add option to save eigenvectors
  void measure_spectral_proj(quad_su3 **conf,theory_pars_t &theory_pars,spectr_proj_meas_pars_t &meas_pars, int iconf,bool conf_created)
  {
    FILE *file=open_file(meas_pars.path,conf_created?"w":"a");
    int neigs = meas_pars.neigs;
    
    complex *charge_cut=nissa_malloc("charge_cut",neigs,complex);
    vector_reset(charge_cut);
    complex *DD_eig_val =nissa_malloc("DD_eig_Val",neigs,complex);
    vector_reset(DD_eig_val);
//    double charge_cut[neigs];
//    double DD_reigs[neigs];

    add_or_rem_stagphases_to_conf(conf); 

    complex **eigvec= nissa_malloc("eigvec",neigs,complex*);
    for(int ieig=0;ieig<neigs;ieig++){
		  eigvec[ieig]=nissa_malloc("eigvec_ieig",(loc_vol+bord_vol)*NCOL,complex);
      vector_reset(eigvec[ieig]);
    }
          
    //loop over hits
    int nhits=meas_pars.nhits;
    for(int hit=0;hit<nhits;hit++)
    {
      verbosity_lv2_master_printf("Evaluating spectral topological charge nhits %d/%d\n",hit+1,nhits);
  
      fill_eigenpart(eigvec,conf,charge_cut,DD_eig_val,meas_pars.neigs,meas_pars.eig_precision);
       
    }
    master_fprintf(file,"%d\t",neigs); 
    for(int ieig=0; ieig<neigs; ++ieig){
      master_fprintf(file,"%.16lg\t%.16lg\t",DD_eig_val[ieig][RE]-0.0001,charge_cut[ieig][RE]);
    }
    master_fprintf(file,"\n");



    add_or_rem_stagphases_to_conf(conf);  
    for(int ieig=0; ieig<neigs; ieig++){
      nissa_free(eigvec[ieig]);
    }
    nissa_free(eigvec);
    nissa_free(charge_cut);
    nissa_free(DD_eig_val);
    close_file(file);

  }
  
  //print pars
  std::string spectr_proj_meas_pars_t::get_str(bool full)
  {
    std::ostringstream os;
    
    os<<"MeasSpectrProj\n";
    os<<base_fermionic_meas_t::get_str(full);
    if(neigs!=def_neigs() or full) os<<" Neigs\t\t=\t"<<neigs<<"\n";
    if(eig_precision!=def_eig_precision() or full) os<<" EigPrecision\t\t=\t"<<eig_precision<<"\n";
    
    return os.str();
  }
}
