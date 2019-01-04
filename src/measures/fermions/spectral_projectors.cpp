#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "dirac_operators/stD/dirac_operator_stD.hpp"
#include "eigenvalues/eigenvalues.hpp"
#include "geometry/geometry_eo.hpp"
#include "geometry/geometry_mix.hpp"
#include "measures/fermions/stag.hpp"

#ifdef USE_THREADS
 #include "routines/thread.hpp"
#endif

#include "spectral_projectors.hpp"

namespace nissa
{
  //This measure will compute the first 'n' eigenvalues (parameter)
  //and eigenvectors of DD^+ in staggered formulation, in order to
  //build an estimate of the topological charge as Q=\sum_i
  //\bar(u)_i gamma5 u_i, where {u_i} are the first n eigenvectors.
  THREADABLE_FUNCTION_7ARG(measure_spectral_proj, complex**,eigvec, quad_su3**,conf,complex*, charge_cut, complex*, DD_eig_val, int,neigs, double,eig_precision, int, wspace_size)
  {
    master_printf("neigs=%d, eig_precision=%.2e, wspace_size=%d\n",neigs,eig_precision,wspace_size);
    
    //wrap the application of DD^+ into an object that can be passed to the eigenfinder
    color *out_tmp_eo[2]={nissa_malloc("out_tmp_EVN",loc_volh+bord_volh,color),nissa_malloc("out_tmp_ODD",loc_volh+bord_volh,color)};
    color *in_tmp_eo[2]={nissa_malloc("in_tmp_EVN",loc_volh+bord_volh,color),nissa_malloc("in_tmp_ODD",loc_volh+bord_volh,color)};
    color *fill_tmp_eo[2]={nissa_malloc("fill_tmp_EVN",loc_volh+bord_volh,color),nissa_malloc("fill_tmp_ODD",loc_volh+bord_volh,color)};
    
    //identity backfield
    quad_u1 *u1b[2]={nissa_malloc("u1b",loc_volh+bord_volh,quad_u1),nissa_malloc("u1b",loc_volh+bord_volh,quad_u1)};
    init_backfield_to_id(u1b);
    
    //Matrix application
    const auto imp_mat=[conf,&in_tmp_eo,&out_tmp_eo](complex *out_lx,complex *in_lx)
      {
	split_lx_vector_into_eo_parts(in_tmp_eo,(color*)in_lx);
	
	evn_apply_stD(out_tmp_eo[EVN],conf,0.00,in_tmp_eo);
	odd_apply_stD(out_tmp_eo[ODD],conf,0.00,in_tmp_eo);
	
	evn_apply_stD_dag(in_tmp_eo[EVN],conf,0.00,out_tmp_eo);
	odd_apply_stD_dag(in_tmp_eo[ODD],conf,0.00,out_tmp_eo);
	
	paste_eo_parts_into_lx_vector((color*)out_lx,in_tmp_eo);
	
	set_borders_invalid(out_lx);
      };
    
    //parameters of the eigensolver
    const bool min_max=0;
    const int mat_size=loc_vol*NCOL;
    const int mat_size_to_allocate=(loc_vol+bord_vol)*NCOL;
    const int niter_max=100000000;
    master_printf("mat_size=%d, mat_size_to_allocate=%d\n",mat_size,mat_size_to_allocate);
    
    //wrap the generation of the test vector into an object that can be passed to the eigenfinder
    const auto filler=[&fill_tmp_eo](complex *out_lx)
      {
      generate_fully_undiluted_eo_source(fill_tmp_eo,RND_GAUSS,-1,0);
      paste_eo_parts_into_lx_vector((color*)out_lx,fill_tmp_eo);
      };
    
    //launch the eigenfinder
    double eig_time=-take_time();
    add_backfield_with_stagphases_to_conf(conf,u1b);
    eigenvalues_of_hermatr_find(eigvec,DD_eig_val,neigs,min_max,mat_size,mat_size_to_allocate,imp_mat,eig_precision,niter_max,filler,wspace_size);
    rem_backfield_with_stagphases_from_conf(conf,u1b);
    
//    complex norm_cut[neigs];
    complex *g5eigvec_i=nissa_malloc("g5eigvec_j",(loc_vol+bord_vol)*NCOL,complex);
    master_printf("\n\nEigenvalues of DD^+:\n");
    for(int ieig=0; ieig<neigs; ++ieig)
      {
      master_printf("lam_%d = (%.16lg,%.16lg)\n",ieig,DD_eig_val[ieig][RE],DD_eig_val[ieig][IM]);
      
      // compute partial sum terms of tr(g5)
      // save 'eigvec[ieig]' in staggered format (i.e., 'fill_tmp_eo'),
      // then multiply it with gamma5 and save the result in
      // 'out_tmp_eo'. The term contributing to the partial sum of tr(g5)
      // will be stored in 'charge_cut[ieig]' as the hermitian product
      // between 'fill_tmp_eo' and 'out_tmp_eo'.
      
      split_lx_vector_into_eo_parts((color**)fill_tmp_eo,(color*)eigvec[ieig]);
      
      //multiply by gamma5
      apply_stag_op(out_tmp_eo,conf,u1b,stag::GAMMA_5,stag::IDENTITY,fill_tmp_eo);
      
      paste_eo_parts_into_lx_vector((color*)g5eigvec_i,out_tmp_eo);
      
      //take hermitian products
      //for(int jeig=ieig; jeig<neigs; ++jeig)
      for(int jeig=ieig; jeig<neigs; ++jeig){
        complex_vector_glb_scalar_prod((double*)&charge_cut[ieig*neigs+jeig],eigvec[jeig],g5eigvec_i,(loc_vol+bord_vol)*NCOL);
        master_printf("u_%d^+ g5 u_%d = (%.16lg,%.16lg)\n",jeig,ieig,charge_cut[ieig*neigs+jeig][RE],charge_cut[ieig*neigs+jeig][IM]);
      }
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
    nissa_free(g5eigvec_i);
  }
  THREADABLE_FUNCTION_END
  
  //measure of spectrally projected components of gamma5 in the staggered formulation
  void measure_spectral_proj(quad_su3 **conf,theory_pars_t &theory_pars,spectr_proj_meas_pars_t &meas_pars,int iconf,bool conf_created)
  {
    /* The format of the file for a single measure is the following:
     *
     * n\t<lam_1>\t<lam_2>\t...\t<lam_n>\t<A_1>\t<A_2>\t...\t<A_n>\t<B_1>\t<B_2>\t...\t<B_n>\n,
     *
     * where:
     *  - 'n' is the number of required eigenvalues (integer)
     *  - 'lam_i' is the i-th eigenvalue of DD^+ (real), 
     *  - A_k = \sumt_{i<=k} u_i^+ g5 u_i (i.e., the k-th partial sum of spectral projections of gamma5; real),
     *  - B_k = \sumt_{i,j<=k} |u_i^+ g5 u_j|^2  (used to estimate the renormalization constant due to 
     *                                            the pseudoscalar current;real and positive)
     */
    FILE *file=open_file(meas_pars.path,conf_created?"w":"a");
    int neigs=meas_pars.neigs;
    
    complex *charge_cut=nissa_malloc("charge_cut",neigs*neigs,complex);
    complex *DD_eig_val=nissa_malloc("DD_eig_Val",neigs,complex);
    vector_reset(charge_cut);
    vector_reset(DD_eig_val);
    
    complex **eigvec=nissa_malloc("eigvec",neigs,complex*);
    for(int ieig=0;ieig<neigs;ieig++)
      {
	eigvec[ieig]=nissa_malloc("eigvec_ieig",(loc_vol+bord_vol)*NCOL,complex);
	vector_reset(eigvec[ieig]);
      }
    
    //loop over hits
    int nhits=meas_pars.nhits;
    for(int hit=0;hit<nhits;hit++)
      {
	verbosity_lv2_master_printf("Evaluating spectral projections of gamma5, nhits %d/%d\n",hit+1,nhits);
	
	measure_spectral_proj(eigvec,conf,charge_cut,DD_eig_val,meas_pars.neigs,meas_pars.eig_precision,meas_pars.wspace_size);
      }
    
    //print the result on file
    master_fprintf(file,"%d\t",neigs);
    for(int ieig=0;ieig<neigs;++ieig)
      master_fprintf(file,"%.16lg\t",DD_eig_val[ieig][RE]);

    double tmp_cum_sum=0.;
    for(int ieig=0;ieig<neigs;++ieig){
      tmp_cum_sum+=charge_cut[ieig*neigs+ieig][RE];
      master_fprintf(file,"%.16lg\t",tmp_cum_sum);
    }
    tmp_cum_sum=0.;
    for(int kcutoff=0; kcutoff<neigs; ++kcutoff){
      tmp_cum_sum+=complex_norm2(charge_cut[kcutoff*neigs+kcutoff]);
      for(int ieig=0; ieig<kcutoff; ++ieig){
        tmp_cum_sum+=2.*complex_norm2(charge_cut[ieig*neigs+kcutoff]);
      }
      master_fprintf(file,"%.16lg\t",tmp_cum_sum);
    }
    master_fprintf(file,"\n");
    
    close_file(file);
    
    for(int ieig=0;ieig<neigs;ieig++)
      nissa_free(eigvec[ieig]);
    nissa_free(eigvec);
    nissa_free(charge_cut);
    nissa_free(DD_eig_val);
  }
  
  //print pars
  std::string spectr_proj_meas_pars_t::get_str(bool full)
  {
    std::ostringstream os;
    
    os<<"MeasSpectrProj\n";
    os<<base_fermionic_meas_t::get_str(full);
    if(neigs!=def_neigs() or full) os<<" Neigs\t\t=\t"<<neigs<<"\n";
    if(eig_precision!=def_eig_precision() or full) os<<" EigPrecision\t\t=\t"<<eig_precision<<"\n";
    if(eig_precision!=def_eig_precision() or full) os<<" WSpaceSize\t\t=\t"<<wspace_size<<"\n";
    
    return os.str();
  }
}
