#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <eigenvalues/eigenvalues.hpp>
#include <dirac_operators/stD/dirac_operator_stD.hpp>

namespace nissa
{
  //computes the spectrum of the overlap operator
  void find_eigenvalues_staggered_DDee(color **eigvec,complex *eigval,int neigs,bool min_max,quad_su3 **conf,color *temp, double mass2, double residue, int wspace_size)
  {
    //Application of the staggered Operator
    const auto imp_mat=[conf,&temp,mass2](complex *out_e,complex *in_e)
      {
	apply_stD2ee_m2((color*)out_e,conf,temp,mass2,(color*)in_e);
      };
    
    const auto filler=[](complex *out_e){ /*currently unimplemented*/};
    
    //parameters of the finder
    const int mat_size=loc_volh*NCOL;
    const int mat_size_to_allocate=(loc_volh+bord_volh)*NCOL;
    const int niter_max=100000000;
    master_printf("mat_size=%d, mat_size_to_allocate=%d\n",mat_size,mat_size_to_allocate);
    
    //precision of the eigenvalues
    double maxerr=sqrt(residue);
    
    verbosity_lv1_master_printf("Starting to search for %d %s eigenvalues of the Staggered operator DD^+ even-projected, with a precision of %lg, and Krylov space size of %d\n",neigs,(min_max?"max":"min"),maxerr,wspace_size);
    
    //find eigenvalues and eigenvectors of the staggered 
    eigenvalues_find((complex**)eigvec,eigval,neigs,min_max,mat_size,mat_size_to_allocate,imp_mat,maxerr,niter_max,filler,wspace_size);
  }
}
