#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif


#include "dirac_operators/overlap/dirac_operator_overlap.hpp"
#include "minmax_eigenvalues.hpp"
#include "eigenvalues/eigenvalues.hpp"
#include "operations/remez/remez_algorithm.hpp"
#include "new_types/rat_approx.hpp"

#ifdef USE_THREADS
 #include "routines/thread.hpp"
#endif

namespace nissa
{
  //measure minmax_eigenvalues
  void measure_minmax_eigenvalues(quad_su3 **conf,theory_pars_t &theory_pars,minmax_eigenvalues_meas_pars_t &meas_pars,int iconf,int conf_created)
  {
    FILE *file=open_file(meas_pars.path,conf_created?"w":"a");
    
    close_file(file);
    
    const double minimum=0.003;
    const double maximum=1.0;
    const int num=-1;
    const int den=2;
    const double tollerance=0.01;
    const double minerr=0.0;
    
    const double t_in=take_time();
    rat_approx_t appr;
    appr.resize(3);
    double res=generate_approx(appr,minimum,maximum,num,den,minerr,tollerance);
    master_printf("Result, res: %lg\n",res);
    appr.master_fprintf_expr(stdout);
    printf("time required= %.10e secs\n",take_time()-t_in);
   
    ////////////////////////////////////////////////////
    // FROM HERE C. BONANNO AND M.CARDINALI OVERLAP ////
    ////////////////////////////////////////////////////

    complex* D_ov_eig_val=nissa_malloc("D_ov_eig_val",meas_pars.neig,complex)
    spincolor** eigvec=nissa_malloc("eigvec",meas_pars.neig,spincolor*);
    for(int ieig=0;ieig<meas_pars.neig;ieig++)
    eigvec[ieig]=nissa_malloc("eig", loc_vol+bord_vol, spincolor);

    

    master_printf("neigs=%d, eig_precision=%.2e\n",meas_pars.neigs,meas_pars.eig_precision);

    //Application of the Overlap Operator
    const auto imp_mat=[conf,M,minerr](complex* out_lx,complex *in_lx)
	{
           apply_overlap((spincolor*)out_lx,conf, theory_pars[0].mass_overlap, minerr, (spincolor*)in_lx);
      	};
    const auto filler=[](complex *out_lx){generate_undiluted_source((spincolor*)out_lx,RND_GAUSS,-1);};
 
    //Parameters of the eigensolver
    const int mat_size=loc_vol*NCOL;
    const int mat_size_to_allocate=(loc_vol+bord_vol)*NCOL;
    const int niter_max=100000000;
    master_printf("mat_size=%d, mat_size_to_allocate=%d\n", mat_size, mat_size_to_allocate);
 
    double eig_time=-take_time();

    //Find eigenvalues and eigenvectors of the overlap
    eigenvalues_of_hermatr_find(eigvec, D_ov_eig_val, meas_pars.neigs, meas_pars.min_max, mat_size, mat_size_to_allocate, imp_mat, meas_pars.eig_precision, niter_max, filler);
 
    master_printf("\n\nEigenvalues of D Overlap:\n");
    for(int ieig=0;ieig<meas_pars.neigs;++ieig)
      {
	 master_printf("%d (%.16lg,%.16lg)\n)",ieig, D_ov_eig_val[RE], D_ov_eig_val[IM]);
      }

    master_printf("\n\n\n");
 
    eig_time+=take_time();
    master_printf("Eigenvalues time: %lg\n", eig_time);
    }
  
  //print
  std::string minmax_eigenvalues_meas_pars_t::get_str(bool full)
  {
    std::ostringstream os;
    
    os<<"MeasMinMaxEigenval\n";
    os<<base_fermionic_meas_t::get_str(full);
    
    return os.str();
  }
}
