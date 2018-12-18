////// METTERE TUTTI GLI INCLUDE ////////

#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "dirac_operators/overlap/dirac_operator_overlap.hpp"
#include "eigenvalues/eigenvalues.hpp"
#include "eig_overlap_operator.hpp"


#ifdef USE_THREADS
 #include "routines/thread.hpp"
#endif

namespace nissa
{ 
  //Written by C. Bonanno and M. Cardinali
  //
  //This measure will compute the first "n" eigenvalues (parameter)
  //and eigenvectors of the Overlap operator D.
  //MAGARI CALCOLIAMO PURE LA CARICA TOPOLOGICA
  //
  THREADABLE_FUNCTION_7ARG(measure_eig_overlap, complex**,eigvec, quad_su3**,conf,complex*, D_ov_eig_val, double, M, double, minerr, int,neigs, double,eig_precision)
 {
  master_printf("neigs=%d, eig_precision=%.2e\n",neigs,eig_precision);

  //Application of the Overlap Operator
  const auto imp_mat=[conf,M,minerr](complex* out_lx,complex *in_lx)
	{
	   apply_overlap((spincolor*)out_lx,conf, M, minerr, (spincolor*)in_lx);
	};
  const auto filler=[](complex *out_lx){generate_undiluted_source((spincolor*)out_lx,RND_GAUSS,-1);};
 
  //Parameters of the eigensolver
  const bool min_max=0;
  const int mat_size=loc_vol*NCOL;
  const int mat_size_to_allocate=(loc_vol+bord_vol)*NCOL;
  const int niter_max=100000000;
  master_printf("mat_size=%d, mat_size_to_allocate=%d\n", mat_size, mat_size_to_allocate);
 
  double eig_time=-take_time();

  //Find eigenvalues and eigenvectors of the overlap
  eigenvalues_of_hermatr_find(eigvec, D_ov_eig_val, neigs, min_max, mat_size, mat_size_to_allocate, imp_mat, eig_precision, niter_max, filler);
 
  master_printf("\n\nEigenvalues of D Overlap:\n");
  for(int ieig=0;ieig<neigs;++ieig)
    {
	 master_printf("%d (%.16lg,%.16lg)\n)",ieig, D_ov_eig_val[RE], D_ov_eig_val[IM]);
    }

  master_printf("\n\n\n");
 
  eig_time+=take_time();
  master_printf("Eigenvalues time: %lg\n", eig_time);
 }
THREADABLE_FUNCTION_END
}
