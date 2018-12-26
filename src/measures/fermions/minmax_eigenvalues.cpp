#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "dirac_operators/overlap/dirac_operator_overlap.hpp"
#include "eigenvalues/eigenvalues.hpp"
#include "geometry/geometry_mix.hpp"
#include "minmax_eigenvalues.hpp"
#include "new_types/rat_approx.hpp"
#include "operations/remez/remez_algorithm.hpp"

#ifdef USE_THREADS
 #include "routines/thread.hpp"
#endif

//////////////////////////////////////////////////////
////      C. BONANNO AND M.CARDINALI OVERLAP      ////
//////////////////////////////////////////////////////

namespace nissa
{
  //measure minmax_eigenvalues
  void measure_minmax_eigenvalues(quad_su3 **conf_eo,theory_pars_t &theory_pars,minmax_eigenvalues_meas_pars_t &meas_pars,int iconf,int conf_created)
  {
    FILE *file=open_file(meas_pars.path,conf_created?"w":"a");
    
    close_file(file);
    
    const double minimum=0.003;
    const double maximum=1.0;
    const int num=-1;
    const int den=2;
    const double tollerance=0.01;
    const double minerr=0.0;
    
    //lx version
    quad_su3 *conf_lx=nissa_malloc("conf_lx",loc_vol+bord_vol,quad_su3);
    paste_eo_parts_into_lx_vector(conf_lx,conf_eo);
    
    const double t_in=take_time();
    rat_approx_t appr;
    appr.resize(3);
    double res=generate_approx(appr,minimum,maximum,num,den,minerr,tollerance);
    master_printf("Result, res: %lg\n",res);
    appr.master_fprintf_expr(stdout);
    printf("time required=%.10e secs\n",take_time()-t_in);
    
    //Parameters of the eigensolver
    const int mat_size=loc_vol*NCOL;
    const int mat_size_to_allocate=(loc_vol+bord_vol)*NCOL;
    const int niter_max=100000000;
    master_printf("mat_size=%d, mat_size_to_allocate=%d\n",mat_size,mat_size_to_allocate);
    
    //allocate
    complex *D_ov_eig_val=nissa_malloc("D_ov_eig_val",meas_pars.neigs,complex);
    spincolor **eigvec=nissa_malloc("eigvec",meas_pars.neigs,spincolor*);
    for(int ieig=0;ieig<meas_pars.neigs;ieig++)
    eigvec[ieig]=nissa_malloc("eig",loc_vol+bord_vol,spincolor);
    
    master_printf("neigs=%d, eig_precision=%.2e\n",meas_pars.neigs,meas_pars.eig_precision);
    
    //consider only the first quark
    int iquark=0;
    if(theory_pars.nflavs()!=1) crash("implemented only for 1 flavor");
    if(theory_pars.quarks[0].discretiz!=ferm_discretiz::OVERLAP) crash("Implemented only for overlap");
    
    //Application of the Overlap Operator
    const auto imp_mat=[conf_lx,&theory_pars,&minerr,iquark](complex* out_lx,complex *in_lx)
	{
	  apply_overlap((spincolor*)out_lx,conf_lx,theory_pars.quarks[iquark].mass_overlap,minerr,(spincolor*)in_lx);
      	};
    const auto filler=[](complex *out_lx){generate_undiluted_source((spincolor*)out_lx,RND_GAUSS,-1);};
    
    double eig_time=-take_time();
    
    //Find eigenvalues and eigenvectors of the overlap
    eigenvalues_find((complex**)eigvec,D_ov_eig_val,meas_pars.neigs,meas_pars.min_max,mat_size,mat_size_to_allocate,imp_mat,meas_pars.eig_precision,niter_max,filler);
    
    master_printf("\n\nEigenvalues of D Overlap:\n");
    for(int ieig=0;ieig<meas_pars.neigs;++ieig)
      master_printf("%d(%.16lg,%.16lg)\n)",ieig,D_ov_eig_val[RE],D_ov_eig_val[IM]);
    
    master_printf("\n\n\n");
    
    eig_time+=take_time();
    master_printf("Eigenvalues time: %lg\n", eig_time);
    
    nissa_free(conf_lx);
  }
  
  //print
  std::string minmax_eigenvalues_meas_pars_t::get_str(bool full)
  {
    std::ostringstream os;
    
    os<<"MeasMinMaxEigenval\n";
    os<<base_fermionic_meas_t::get_str(full);
    if(neigs!=def_neigs() or full) os<<" Neigs\t\t=\t"<<neigs<<"\n";
    if(eig_precision!=def_eig_precision() or full) os<<" EigPrecision\t\t=\t"<<eig_precision<<"\n";
    if(wspace_size!=def_wspace_size() or full) os<<" WSpaceSize\t\t=\t"<<wspace_size<<"\n";
    if(min_max!=def_min_max() or full) os<<" MinMax\t\t=\t"<<min_max<<"\n";
    
    return os.str();
  }
}
