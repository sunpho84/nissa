#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "linalgs/reduce.hpp"
#include "measures/fermions/tm_corr_op.hpp"

namespace nissa
{
  void tm_corr_op::undiluted_meson_contr(complex* contr,
					 spincolor *bw,
					 spincolor *fw,
					 const int& igamma,
					 const int source_coord)
  {
    /// Local storage
    complex *loc_contr=get_reducing_buffer<complex>(locVol);
    
    /// Gamma matrix
    dirac_matr g=base_gamma[igamma]*base_gamma[5];
    
    NISSA_PARALLEL_LOOP(ivol,0,locVol)
      {
	spincolor temp;
	unsafe_dirac_prod_spincolor(temp,&g,fw[ivol]);
	spincolor_scalar_prod(loc_contr[ivol],bw[ivol],temp);
      }
    NISSA_PARALLEL_LOOP_END;
    THREAD_BARRIER();
    
    /// Temporary contraction
    complex unshiftedGlbContr[glbSize[0]];
    glb_reduce(unshiftedGlbContr,loc_contr,locVol,glbSize[0],locSize[0],glbCoordOfLoclx[0][0]);
    
    for(int glb_t=0;glb_t<glbSize[0];glb_t++)
      {
	/// Distance from source
	const int dt=
	  (glb_t-source_coord+glbSize[0])%glbSize[0];
	
	complex_copy(contr[dt],unshiftedGlbContr[glb_t]);
      }
    
  }
}
