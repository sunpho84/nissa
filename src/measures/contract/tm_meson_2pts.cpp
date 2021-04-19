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
					 const GlbCoord& source_coord)
  {
    /// Local storage
    complex *loc_contr=get_reducing_buffer<complex>(locVol());
    
    /// Gamma matrix
    dirac_matr g=base_gamma[igamma]*base_gamma[5];
    
    NISSA_PARALLEL_LOOP(ivol,0,locVol)
      {
	spincolor temp;
	unsafe_dirac_prod_spincolor(temp,&g,fw[ivol.nastyConvert()]);
	spincolor_scalar_prod(loc_contr[ivol.nastyConvert()],bw[ivol.nastyConvert()],temp);
      }
    NISSA_PARALLEL_LOOP_END;
    THREAD_BARRIER();
    
    /// Temporary contraction
    complex unshiftedGlbContr[glbTimeSize.nastyConvert()];
    glb_reduce(unshiftedGlbContr,loc_contr,locVol(),glbTimeSize(),locTimeSize(),glbCoordOfLoclx(LocLxSite(0),timeDirection)());
    
    for(GlbCoord glb_t=0;glb_t<glbTimeSize;glb_t++)
      {
	/// Distance from source
	const GlbCoord dt=
	  (glb_t-source_coord+glbTimeSize)%glbTimeSize;
	
	complex_copy(contr[dt()],unshiftedGlbContr[glb_t.nastyConvert()]);
      }
    
  }
}
