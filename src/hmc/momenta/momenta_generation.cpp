#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "base/random.hpp"
#include "base/vectors.hpp"
#include "dirac_operators/momenta/MFACC.hpp"
#include "geometry/geometry_eo.hpp"
#include "hmc/gauge/MFACC_fields.hpp"
#include "inverters/momenta/cg_invert_MFACC.hpp"
#include "inverters/momenta/cgm_invert_MFACC.hpp"
#include "linalgs/linalgs.hpp"
#include "new_types/rat_approx.hpp"
#include "new_types/su3_op.hpp"
#include "operations/remez/remez_algorithm.hpp"
#include "routines/ios.hpp"
#include "threads/threads.hpp"

namespace nissa
{
  //generate momenta using guassian hermitean matrix generator
  void generate_hmc_momenta(eo_ptr<quad_su3> H)
  {
    for(int par=0;par<2;par++)
      {
	NISSA_PARALLEL_LOOP(ieo,0,locVolh)
	  for(int mu=0;mu<NDIM;mu++)
	    herm_put_to_gauss(H[par][ieo][mu],&(loc_rnd_gen[loclx_of_loceo[par][ieo]]),1);
	NISSA_PARALLEL_LOOP_END;
	
	set_borders_invalid(H[par]);
      }
  }
  //similar for lx
  void generate_hmc_momenta(quad_su3* H)
  {
    
    NISSA_PARALLEL_LOOP(ivol,0,locVol)
      for(int mu=0;mu<NDIM;mu++)
	herm_put_to_gauss(H[ivol][mu],&(loc_rnd_gen[ivol]),1);
    NISSA_PARALLEL_LOOP_END;
    
    set_borders_invalid(H);
  }
  
  //generate momenta using guassian hermitian matrix generator
  void generate_hmc_momenta_with_FACC(quad_su3* H,quad_su3* conf,rat_approx_t* rat_exp_H,double kappa,double residue)
  {
    
    //temporary for inversion
    su3 *in=nissa_malloc("in",locVol+bord_vol,su3);
    su3 *out=nissa_malloc("out",locVol+bord_vol,su3);
    su3 *tmp=nissa_malloc("tmp",locVol+bord_vol,su3);
    
    for(int mu=0;mu<NDIM;mu++)
      {
	//fill the vector randomly
	NISSA_PARALLEL_LOOP(ivol,0,locVol)
	  herm_put_to_gauss(in[ivol],&(loc_rnd_gen[ivol]),1);
	NISSA_PARALLEL_LOOP_END;
	set_borders_invalid(in);
	
	//compute the norm
	double norm;
	double_vector_glb_scalar_prod(&norm,(double*)in,(double*)in,locVol*sizeof(su3)/sizeof(double));
	
	//invert
	summ_src_and_all_inv_MFACC_cgm(out,conf,kappa,rat_exp_H,1000000,residue,in);
	
	//try to compute the norm*D
        inv_MFACC_cg(tmp,NULL,conf,kappa,10000000,residue,out);
	double norm_reco;
	double_vector_glb_scalar_prod(&norm_reco,(double*)out,(double*)tmp,locVol*sizeof(su3)/sizeof(double));
	master_printf("Norm: %16.16lg, norm_reco: %16.16lg, relative error: %lg\n",sqrt(norm),sqrt(norm_reco),sqrt(norm/norm_reco)-1);
	
	//store the vector
	NISSA_PARALLEL_LOOP(ivol,0,locVol)
	  su3_copy(H[ivol][mu],out[ivol]);
	NISSA_PARALLEL_LOOP_END;
	set_borders_invalid(H);
      }
    
    nissa_free(in);
    nissa_free(out);
    nissa_free(tmp);
  }
  
  //generate momenta needed for Fourier acceleration
  void generate_MFACC_momenta(su3** pi,int naux_fields,quad_su3* conf,rat_approx_t* rat_exp_H,double kappa,double residue)
  {
    verbosity_lv1_master_printf("Generating Fourier acceleration momenta\n");
    
    //allocate gaussian field
    su3 *V=nissa_malloc("V",locVol+bord_vol,su3);
    
    double act=0;
    for(int id=0;id<naux_fields;id++)
      {
        //generate gaussianly
        generate_MFACC_field(pi[id]);
	
	//compute act
	double temp;
	double_vector_glb_scalar_prod(&temp,(double*)(pi[id]),(double*)(pi[id]),locVol*sizeof(su3)/sizeof(double));
	act+=temp/2;
	
	//multiply by sqrt
	summ_src_and_all_inv_MFACC_cgm(V,conf,kappa,rat_exp_H,1000000,residue,pi[id]);
	inv_MFACC_cg(pi[id],NULL,conf,kappa,10000000,residue,V);
      }
    verbosity_lv1_master_printf("Act internal: %+16.16lg\n",act);
    
    nissa_free(V);
  }
}
