#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "geometry/geometry_lx.hpp"
#include "linalgs/reduce.hpp"
#include "new_types/complex.hpp"
#include "new_types/su3.hpp"

#include "measures/contract/baryon.hpp"

namespace nissa
{
  /// Compute G*gammaX
  dirac_matr get_CgX(const int X)
  {
    dirac_matr CgX;
    
    dirac_matr g2g4,C;
    dirac_prod(&g2g4,base_gamma+2,base_gamma+4);
    dirac_prod_idouble(&C,&g2g4,1);
    dirac_prod(&CgX,&C,base_gamma+X);
    
    return CgX;
  }
  
  void compute_baryon_2pts_proj_contr(complex* contr,
				      const int& igSo,
				      const int& igSi,
				      spincolor** Q1,
				      spincolor** Q2,
				      spincolor** Q3,
				      const int source_coord,
				      const double& temporal_bc)
  {
    /// Pack the three propagators
    spincolor** Q[3]={Q1,Q2,Q3};
    
    /// Index of the permutation
    const int eps_sign[2]={1,-1};
    
    /// List the two other terms for each color index
    const std::array<std::array<int,2>,3> eps_id={{{1,2}, {2,0}, {0,1}}};
    
    /// Number of different contributions Id and gamma0
    const int nIdg0=2;
    
    /// Number of wicks contractions
    const int nWicks=2;
    
    //When taking the adjoint interpolating operator, we must
    //include the sign of g0 Cg^\dagger g0.
    const auto& g=base_gamma;
    
    const dirac_matr Cg_so=herm(g[4]*get_CgX(igSo)*g[4]);
    const dirac_matr Cg_si=get_CgX(igSi);
    
    //Precompute the factor to be added
    spinspin fact;
    for(int sp_si=0;sp_si<NDIRAC;sp_si++)
      for(int sp_so=0;sp_so<NDIRAC;sp_so++)
	{
	  complex& f=fact[sp_si][sp_so];
	  unsafe_complex_prod(f,Cg_si.entr[sp_si],Cg_so.entr[sp_so]);
	}
    
    //allocate loc storage
    complex *loc_contr=get_reducing_buffer<complex>(loc_vol*nIdg0*nWicks);
    vector_reset(loc_contr);
    
    NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
      {
	/// Distance from source
	const int dt=(glb_coord_of_loclx[ivol][0]-source_coord+glb_size[0])%glb_size[0];
	
	/// Determine whether we are in the first half
	const bool first_half=(dt<=glb_size[0]/2);
	
	//Compute the projector, gi*gj*(1 or g0)
	dirac_matr proj[nIdg0];
	const int g_of_id_g0[nIdg0]={0,4};
	for(int idg0=0;idg0<nIdg0;idg0++)
	  proj[idg0]=g[igSi]*g[igSo]*g[g_of_id_g0[idg0]];
	
	//Takes a slice
	su3spinspin p1,p2,p3;
	su3spinspin* p[3]={&p1,&p2,&p3};
	for(int i=0;i<3;i++)
	  for(int sp_so=0;sp_so<NDIRAC;sp_so++)
	    for(int sp_si=0;sp_si<NDIRAC;sp_si++)
	      for(int co_so=0;co_so<NCOL;co_so++)
		for(int co_si=0;co_si<NCOL;co_si++)
		  complex_copy((*p)[i][co_si][co_so][sp_si][sp_so],Q[i][co_so+NCOL*sp_so][ivol][sp_si][co_si]);
	
	//Color source
	for(int b_so=0;b_so<NCOL;b_so++)
	  //Color sink
	  for(int b_si=0;b_si<NCOL;b_si++)
	    {
	      //Dirac source
	      for(int al_so=0;al_so<NDIRAC;al_so++)
		for(int be_so=Cg_so.pos[al_so],
		      //Dirac sink
		      al_si=0;al_si<NDIRAC;al_si++)
		  {
		    const int be_si=Cg_si.pos[al_si];
		    
		    complex AC_proj[2]={};
		    
		    for(int iperm_so=0;iperm_so<2;iperm_so++)
		      for(int iperm_si=0;iperm_si<2;iperm_si++)
			{
			  const int c_so=eps_id[b_so][1-iperm_so],a_so=eps_id[b_so][iperm_so];
			  const int c_si=eps_id[b_si][1-iperm_si],a_si=eps_id[b_si][iperm_si];
			  
			  for(int ga_so=0;ga_so<NDIRAC;ga_so++)
			    for(int idg0=0;idg0<nIdg0;idg0++)
			      {
				const int ga_si=proj[idg0].pos[ga_so];
				
				//In practice, only in the second half and for the g0 we have a minus
				
				const int sign_idg0=(first_half or idg0==0)?+1:-1;
				const int sign_tot=eps_sign[iperm_so]*eps_sign[iperm_si]*sign_idg0;
				
				const int sp1_si[2]={al_si,ga_si};
				const int sp3_si[2]={ga_si,al_si};
				
				const int co1_si[2]={a_si,c_si};
				const int co3_si[2]={c_si,a_si};
				
				for(int iWick=0;iWick<nWicks;iWick++)
				  {
				    complex AC;
				    unsafe_complex_prod(AC,p1[co1_si[iWick]][a_so][sp1_si[iWick]][al_so],p3[co3_si[iWick]][c_so][sp3_si[iWick]][ga_so]);
				    complex_prodassign_double(AC,sign_tot);
				    complex_summ_the_prod(AC_proj[iWick],AC,proj[idg0].entr[ga_so]);
				  }
			      }
			}
		    
		    /// Local time
		    const int loc_t=loc_coord_of_loclx[ivol][0];
		    
		    /// Local spatial id
		    const int loc_ispat=ivol%loc_spat_vol;
		    
		    for(int iWick=0;iWick<nWicks;iWick++)
		      {
			/// Index in the local elements: we put space most internally, and local time most externally
			const int iloc=loc_ispat+loc_spat_vol*(iWick+nWicks*loc_t);
			
			complex term;
			unsafe_complex_prod(term,AC_proj[iWick],fact[al_si][al_so]);
			complex_summ_the_prod(loc_contr[iloc],term,p2[b_si][b_so][be_si][be_so]);
		      }
		  }
	    }
      }
    NISSA_PARALLEL_LOOP_END;
    THREAD_BARRIER();
    
    /// Number of local elements
    const int nloc=nWicks*loc_vol;
    
    /// Number of slices to be taken
    const int nslices=nWicks*glb_size[0];
    
    /// Number of local slices referring to the elements
    const int nloc_slices=nWicks*loc_size[0];
    
    /// Offset of each local slice w.r.t. globals
    const int loc_offset=nWicks*glb_coord_of_loclx[0][0];
    
    complex unshifted_glb_contr[glb_size[0]*nWicks];
    glb_reduce(unshifted_glb_contr,loc_contr,nloc,nslices,nloc_slices,loc_offset);
    
    for(int glb_t=0;glb_t<glb_size[0];glb_t++)
      for(int iWick=0;iWick<nWicks;iWick++)
	{
	  /// Distance from source
	  const int dt=
	    (glb_t-source_coord+glb_size[0])%glb_size[0];
	  
	  /// Input index
	  const int iin=iWick+nWicks*glb_t;
	  
	  /// Argument of the phase
	  const double arg=3*temporal_bc*M_PI*dt/glb_size[0];
	  
	  /// Phase
	  const complex phase={cos(arg),sin(arg)};
	  
	  complex_summ_the_prod(contr[iWick+nWicks*dt],unshifted_glb_contr[iin],phase);
	}
  }
  
  void compute_nucleon_2pts_contr(complex* contr,
				  spincolor** Ql,
				  spincolor** Qd,
				  const int source_coord,
				  const double& temporal_bc)
  {
    /// Uncombined contraction, with the two Wick contractions separately
    complex contr_per_Wick[glb_size[0]*2];
    compute_baryon_2pts_proj_contr(contr_per_Wick,5,5,Ql,Qd,Ql,source_coord,temporal_bc);
    
    for(int t=0;t<glb_size[0];t++)
      complex_subt(contr[t],contr_per_Wick[0+2*t],contr_per_Wick[1+2*t]);
  }
}
