#ifdef HAVE_CONFIG_H
# include "config.hpp"
#endif

#include <dirac_operators/overlap/dirac_operator_overlap.hpp>
#include <eigenvalues/eigenvalues.hpp>
#include <random/randomGenerate.hpp>

namespace nissa
{
  //computes the spectrum of the overlap operator
  void find_eigenvalues_overlap(spincolor **eigvec,complex *eigval,int neigs,bool min_max,quad_su3 *conf,rat_approx_t& appr,double residue,double mass_overlap,double mass,int wspace_size)
  {
    //Application of the overlap Operator
    const auto imp_mat=[conf,&appr,residue,mass_overlap,mass](complex *out_lx,complex *in_lx)
      {
	apply_overlap((spincolor*)out_lx,conf,&appr,residue,mass_overlap,mass,(spincolor*)in_lx);
      };
    
    const auto filler=[](complex *out_lx)
      {
	generate_undiluted_source((spincolor*)out_lx,RND_GAUSS,-1);
      };
    
    //parameters of the finder
    const int mat_size=locVol.nastyConvert()*NCOL*NDIRAC;
    const int mat_size_to_allocate=locVolWithBord.nastyConvert()*NCOL*NDIRAC;
    const int niter_max=100000000;
    master_printf("mat_size=%d, mat_size_to_allocate=%d\n",mat_size,mat_size_to_allocate);
    
    //precision of the eigenvalues
    double maxerr=sqrt(residue);
    
    verbosity_lv1_master_printf("Starting to search for %d %s eigenvalues of the overlap operator, with a precision of %lg, and Krylov space size of %d\n",neigs,(min_max?"max":"min"),maxerr,wspace_size);
    
    //find eigenvalues and eigenvectors of the overlap
    eigenvalues_find((complex**)eigvec,eigval,neigs,min_max,mat_size,mat_size_to_allocate,imp_mat,maxerr,niter_max,filler,wspace_size);
  }
}
