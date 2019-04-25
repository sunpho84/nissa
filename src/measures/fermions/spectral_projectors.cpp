#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "dirac_operators/stD/dirac_operator_stD.hpp"
#include "eigenvalues/eigenvalues.hpp"
#include "eigenvalues/eigenvalues_staggered.hpp"
#include "geometry/geometry_eo.hpp"
#include "geometry/geometry_mix.hpp"
#include "measures/fermions/stag.hpp"

#ifdef USE_THREADS
 #include "routines/thread.hpp"
#endif

#include "spectral_projectors.hpp"

namespace nissa
{

  // This measure will compute the first 'n' eigenvalues (parameter)
  // and eigenvectors of the iD operator in staggered formulation, in order to
  // build an estimate of the topological susceptibility.
  // refs:  https://arxiv.org/pdf/1008.0732.pdf for the susceptibility formula,
  //        https://arxiv.org/pdf/0912.2850.pdf for the 2^(d/2) overcounting.
  THREADABLE_FUNCTION_7ARG(measure_iD_spectrum, color**,eigvec, quad_su3**,conf,complex*, charge_cut, complex*, eigval, int,neigs, double,eig_precision, int,wspace_size)
  {
    //parameters of the eigensolver
    const bool min_max=0;
    
    //identity backfield
    quad_u1 *u1b[2]={nissa_malloc("u1b",loc_volh+bord_volh,quad_u1),nissa_malloc("u1b",loc_volh+bord_volh,quad_u1)};
    init_backfield_to_id(u1b);

    //temporary vectors
    color *tmpvec_eo[2]={nissa_malloc("tmpvec_eo_EVN",loc_volh+bord_volh,color),nissa_malloc("tmpvec_eo_ODD",loc_volh+bord_volh,color)};

    //results of the g5 application
    color *eigvec_g5_eo[2]={nissa_malloc("eigvec_g5_EVN",loc_volh+bord_volh,color),nissa_malloc("eigvec_g5_ODD",loc_volh+bord_volh,color)};
    color *eigvec_g5_lx=nissa_malloc("eigvec_g5",loc_vol+bord_vol,color);
    
    //launch the eigenfinder
    double eig_time=-take_time();
    find_eigenvalues_staggered_iD(eigvec,eigval,neigs,min_max,conf,u1b,eig_precision,wspace_size);
    
    verbosity_lv1_master_printf("\n\nEigenvalues of staggered iD operator:\n");
    for(int ieig=0;ieig<neigs;ieig++)
    {
      master_printf("lam_%d = (%.16lg,%.16lg)\n",ieig,eigval[ieig][RE],eigval[ieig][IM]);

      // compute terms u_j^+ g5 u_i
      // convert 'eigvec[ieig]' in staggered format ('in_tmp_eo'),
      // then multiply it with gamma5 and save the result in
      // 'out_tmp_eo'. The term corresponding to u_j^+ g5 u_i
      // will be stored in 'charge_cut[ieig*neigs+jeigs]' as the hermitian product
      // between 'eigvec[jeig]' and 'out_tmp_eo[EVN]'.

      //multiply by gamma5
      split_lx_vector_into_eo_parts(tmpvec_eo,eigvec[ieig]);
      apply_stag_op(eigvec_g5_eo,conf,u1b,stag::GAMMA_5,stag::IDENTITY,tmpvec_eo);
      paste_eo_parts_into_lx_vector(eigvec_g5_lx,eigvec_g5_eo);

      //take hermitian products 
      for(int jeig=ieig;jeig<neigs;jeig++)
    {
      complex_vector_glb_scalar_prod(charge_cut[ieig*neigs+jeig],(complex*)eigvec[jeig],(complex*)eigvec_g5_lx,loc_vol*sizeof(color)/sizeof(complex));
      verbosity_lv1_master_printf("u_%d^+ g5 u_%d = (%.16lg,%.16lg)\n",jeig,ieig,charge_cut[ieig*neigs+jeig][RE],charge_cut[ieig*neigs+jeig][IM]);
    }

    }
    verbosity_lv2_master_printf("\n\n\n");
    
    eig_time+=take_time();
    verbosity_lv1_master_printf("Eigenvalues time: %lg\n",eig_time);
    
    nissa_free(tmpvec_eo[EVN]);
    nissa_free(tmpvec_eo[ODD]);
    nissa_free(eigvec_g5_eo[EVN]);
    nissa_free(eigvec_g5_eo[ODD]);
    nissa_free(eigvec_g5_lx);
    nissa_free(u1b[0]);
    nissa_free(u1b[1]);
  }
  THREADABLE_FUNCTION_END


  //measure of spectrally projected components of gamma5 in the staggered formulation using iD
  void measure_spectral_proj(quad_su3 **conf,theory_pars_t &theory_pars,spectr_proj_meas_pars_t &meas_pars,int iconf,bool conf_created)
  {
    /* The format of the file for a single measure is the following:
     *
     * n\t<lam_1>\t<lam_2>\t...\t<lam_n>\t<A_1>\t<A_2>\t...\t<A_n>\t<B_1>\t<B_2>\t...\t<B_n>\n,
     *
     * where:
     *  - 'n' is the number of required eigenvalues (integer)
     *  - 'lam_i' is the i-th eigenvalue of the Adams operator (real, ordered from small to large in magnitude),
     *  - A_k = \sumt_{i<=k} u_i^+ g5 u_i (i.e., the k-th partial sum of spectral projections of gamma5; real),
     *  - B_k = \sumt_{i,j<=k} |u_i^+ g5 u_j|^2  (used to estimate the renormalization constant due to
     *                                            the pseudoscalar current;real and positive)
     */
    int neigs=meas_pars.neigs;
   
    // allocate auxiliary vectors 
		quad_su3 *conf_lx=nissa_malloc("conf_lx",loc_vol+bord_vol,quad_su3);
		paste_eo_parts_into_lx_vector(conf_lx,conf);
		quad_su3 *conf_eo[2]={nissa_malloc("conf_eo_EVN",loc_volh+bord_volh,quad_su3),nissa_malloc("conf_eo_ODD",loc_volh+bord_volh,quad_su3)};
		split_lx_vector_into_eo_parts(conf_eo,conf_lx);
    complex *charge_cut=nissa_malloc("charge_cut",neigs*neigs,complex);
    complex *eigval=nissa_malloc("DD_eig_Val",neigs,complex);
    double *cum_sumA=nissa_malloc("cum_sumA",meas_pars.neigs+1,double);
    double *cum_sumB=nissa_malloc("cum_sumB",meas_pars.neigs+1,double);
    color **eigvec=nissa_malloc("eigvec",neigs,color*);
    for(int ieig=0;ieig<neigs;ieig++)
      eigvec[ieig]=nissa_malloc("eigvec_ieig",loc_vol+bord_vol,color);
   

    //loop on smooth
    int nsmooth=0;
    bool finished=true;
		do
			{
    verbosity_lv1_master_printf("Measuring spectral projectors for nsmooth %d/%d\n",nsmooth,meas_pars.smooth_pars.nsmooth());
    // reset vectors
    vector_reset(charge_cut);
    vector_reset(eigval);
    vector_reset(cum_sumA);
    vector_reset(cum_sumB);
    for(int ieig=0;ieig<neigs;ieig++)
      vector_reset(eigvec[ieig]);

    measure_iD_spectrum(eigvec,conf_eo,charge_cut,eigval,meas_pars.neigs,meas_pars.eig_precision,meas_pars.wspace_size);
    
    //print the result on file
    verbosity_lv2_master_printf("\n\nPartial sums for spectral projectors:\n\nk\t\t\teig\t\t\tA_k\t\t\tB_k\n");
    // vectors storing A_k and B_k partial sums, offset by 1 for convenience 
    
    FILE *file=open_file(meas_pars.path,conf_created?"w":"a");

    master_fprintf(file,"%d\t%d\t%d\t",iconf,nsmooth,neigs);
    for(int ieig=0;ieig<neigs;++ieig)
      master_fprintf(file,"%.16lg\t",eigval[ieig][RE]);
    
    for(int ieig=0;ieig<neigs;ieig++)
      {
  cum_sumA[1+ieig]=cum_sumA[ieig]+charge_cut[ieig*neigs+ieig][RE];
  master_fprintf(file,"%.16lg\t",cum_sumA[1+ieig]);
      }

    for(int kcutoff=0;kcutoff<neigs;kcutoff++)
      {
      cum_sumB[1+kcutoff]=cum_sumB[kcutoff]+complex_norm2(charge_cut[kcutoff*neigs+kcutoff]); //diagonal part
      for(int ieig=0;ieig<kcutoff;ieig++)
        cum_sumB[1+kcutoff]+=2.0*complex_norm2(charge_cut[ieig*neigs+kcutoff]); //offdiagonal part

      master_fprintf(file,"%.16lg\t",cum_sumB[1+kcutoff]);
      }
    master_fprintf(file,"\n");
    
    close_file(file);

    for(int ieig=0; ieig<neigs; ieig++){
      verbosity_lv2_master_printf("%d\t%.16lg\t%.16lg\t%.16lg\n",ieig,eigval[ieig][RE],cum_sumA[1+ieig],cum_sumB[1+ieig]);
    }
    verbosity_lv2_master_printf("\n\n");

		//proceeds with smoothing
		if(meas_pars.smooth_pars.next_nsmooth_meas(nsmooth) > meas_pars.smooth_pars.nsmooth())
			break;

		finished=smooth_lx_conf_until_next_meas(conf_lx,meas_pars.smooth_pars,nsmooth);
		split_lx_vector_into_eo_parts(conf_eo,conf_lx);
      }
    while(not finished);
   
    // deallocating vectors 
    for(int ieig=0;ieig<neigs;ieig++)
      nissa_free(eigvec[ieig]);

    nissa_free(conf_lx);
    nissa_free(conf_eo[EVN]);
    nissa_free(conf_eo[ODD]);
    nissa_free(eigvec);
    nissa_free(charge_cut);
    nissa_free(eigval);
    nissa_free(cum_sumA);
    nissa_free(cum_sumB);
  }
  
  //print pars
  std::string spectr_proj_meas_pars_t::get_str(bool full)
  {
    std::ostringstream os;
    
    os<<"MeasSpectrProj\n";
    os<<base_fermionic_meas_t::get_str(full);
    if(neigs!=def_neigs() or full) os<<" Neigs\t\t=\t"<<neigs<<"\n";
    if(eig_precision!=def_eig_precision() or full) os<<" EigPrecision\t\t=\t"<<eig_precision<<"\n";
    if(wspace_size!=def_wspace_size() or full) os<<" WSpaceSize\t\t=\t"<<wspace_size<<"\n";
    os<<smooth_pars.get_str(full);
    
    return os.str();
  }
}
