#ifdef HAVE_CONFIG_H
# include <config.hpp>
#endif

#include <base/field.hpp>
#include "base/vectors.hpp"
#include "communicate/borders.hpp"
#include "communicate/edges.hpp"
#include "geometry/geometry_lx.hpp"
#include "new_types/su3_op.hpp"
#include "routines/ios.hpp"
#include "threads/threads.hpp"

#include <string.h>

//implementation of hep-lat/0607006, Capitani et al. Appendix A.4
//normalization taken from 1011.2711
namespace nissa
{
  template <typename...Tail>
  void sme(LxField<su3>& lOut,
	   const double& alpha,
	   const std::function<LxField<su3>&(int,int,int,int)>& V, // needs to pass through std::function or nvcc complains badly
	   const LxField<su3>& U,
	   const int& mu,
	   const Tail&...a)
  {
    lOut.reset();
    
    for(int sigma=0;sigma<NDIM;sigma++)
      {
	bool toDo=true;
	for(const auto& i : {mu,a...})
	  toDo&=(i!=sigma);
	
	if(toDo)
	  PAR(0,locVol,
	      CAPTURE(mu,
		      sigma,
		      A=V(sigma,mu,a...).getReadable(),
		      B=V(mu,a...,sigma).getReadable(),
		      TO_WRITE(lOut)),
	      iVol,
	      {
		su3 temp;
		
		const int muUp=loclxNeighup[iVol][mu];
		
		unsafe_su3_prod_su3(temp,A[iVol],B[loclxNeighup[iVol][sigma]]);
		su3_summ_the_prod_su3_dag(lOut[iVol],temp,A[muUp]);
		
		const int sigmaDw=loclxNeighdw[iVol][sigma];
		const int sigmaDwMuUp=loclxNeighdw[muUp][sigma];
		
		unsafe_su3_dag_prod_su3(temp,A[sigmaDw],B[sigmaDw]);
		su3_summ_the_prod_su3(lOut[iVol],temp,A[sigmaDwMuUp]);
	      });
      }
    
    PAR(0,
	locVol,
	CAPTURE(TO_WRITE(lOut),
		TO_READ(U),
		alpha),
	iVol,
      {
	su3 temp1,temp2;
	
	unsafe_su3_prod_su3_dag(temp1,lOut[iVol],U[iVol]);
	unsafe_su3_traceless_anti_hermitian_part(temp2,temp1);
	
	su3 Q;
	unsafe_su3_traceless_anti_hermitian_part(Q,temp2);
	su3_prodassign_idouble(Q,-alpha);
	
	safe_hermitian_exact_i_exponentiate(Q,Q);
	unsafe_su3_prod_su3(lOut[iVol],Q,U[iVol]);
      });
  }
  
  /// Hex smear the conf
  void hex_smear_conf(LxField<quad_su3>& sm_conf,
		      const LxField<quad_su3>& conf,
		      const double& alpha1,
		      const double& alpha2,
		      const double& alpha3)
  {
    const int nDecLevels=4;
    const int nDecsPerLevel[]=
      {NDIM,
       NDIM*(NDIM-1),
       NDIM*(NDIM-1),
       NDIM};
    
    const int nDecTot=
		nDecsPerLevel[0]+nDecsPerLevel[1]+nDecsPerLevel[2]+nDecsPerLevel[3];
    std::vector<LxField<su3>> linksAllDecs(nDecTot,{"Dec",WITH_HALO_EDGES});
    
    std::vector<LxField<su3>*> linksDecsPerLev(nDecLevels);
    linksDecsPerLev[0]=&linksAllDecs[0];
    for(int iDecLev=1;iDecLev<nDecLevels;iDecLev++)
      linksDecsPerLev[iDecLev]=linksDecsPerLev[iDecLev-1]+nDecsPerLevel[iDecLev-1];
    
    auto communicateDecsForLev=
      [linksDecsPerLev,
       nDecsPerLevel](const int& iLev)
      {
	for(int iDec=0;iDec<nDecsPerLevel[iLev];iDec++)
	  {
	    LxField<su3>& o=linksDecsPerLev[iLev][iDec];
	    o.invalidateHalo();
	    o.updateEdges();
	  }
      };
    
    /////////////////////////////////////////////////////////////////
    
    /// Level0
    
    const auto links0=
      [links0List=linksDecsPerLev[0]](const int& mu,
				      const int& =0,
				      const int& =0,
				      const int& =0) -> LxField<su3>& //swallow all other dirs
      {
	return
	  *links0List[mu];
      };
    
    for(int mu=0;mu<NDIM;mu++)
      PAR(0,
	  locVol,
	  CAPTURE(mu,
		  TO_READ(conf),
		  u=links0(mu).getWritable()),
	  A,
	  {
	    su3_copy(u[A],conf[A][mu]);
	  });
    
    communicateDecsForLev(0);
    
    /////////////////////////////////////////////////////////////////
    
    /// Level1
    
    const auto links1=
      [links1List=linksDecsPerLev[1]](const int& mu,
				      const int& nu,
				      const int& rho,
				      const int& =0) -> LxField<su3>&
      {
	const int e=15^((1<<mu)+(1<<nu)+(1<<rho));
	const int sigma=(e==2)+2*(e==4)+3*(e==8);
	const int iSigma=sigma-(sigma>mu); /// Subtract 1 if sigma>mu since we don't store mu==sigma
	
	return
	  *links1List[iSigma+(NDIM-1)*mu];
      };
    
    for(int mu=0;mu<NDIM;mu++)
      for(int iNu=0;iNu<NDIM-1;iNu++)
	for(int iRho=0;iRho<NDIM-2;iRho++)
	  {
	    const int nu=perpDirs[mu][iNu];
	    const int rho=perp2Dirs[mu][iNu][iRho];
	    
	    if(rho>nu)
	      sme(links1(mu,nu,rho),alpha3/2,
		  links0,
		  links0(mu),
		  mu,nu,rho);
	  }
    
    communicateDecsForLev(1);
    
    /////////////////////////////////////////////////////////////////
    
    //Level2
    
    const auto links2=
      [links2List=linksDecsPerLev[2]](const int& mu,
				      const int& nu,
				      const int& =0,
				      const int& =0) -> LxField<su3>&
      {
	const int iNu=
	  nu-(nu>mu); /// Subtract 1 if sigma>mu since we don't store mu==sigma
	
	return
	  *links2List[iNu+(NDIM-1)*mu];
      };
    
    for(int mu=0;mu<NDIM;mu++)
      for(int iNu=0;iNu<NDIM-1;iNu++)
	{
	  const int nu=perpDirs[mu][iNu];
	  sme(links2(mu,nu),alpha2/4,links1,links0(mu),mu,nu,0);
	}
    
    communicateDecsForLev(2);
    
    /////////////////////////////////////////////////////////////////
    
    /// Level3
    
    const auto links3=
      [links3List=linksDecsPerLev[3]](const int& mu,
				      const int& =0,
				      const int& =0,
				      const int& =0) -> LxField<su3>&
      {
	return
	  *links3List[mu];
      };
    
    for(int mu=0;mu<NDIM;mu++)
      sme(links3(mu),alpha1/6,links2,links0(mu),mu,0,0);
    
    /////////////////////////////////////////////////////////////////
    
    for(int mu=0;mu<NDIM;mu++)
      PAR(0,
	  locVol,
	  CAPTURE(mu,
		  TO_WRITE(sm_conf),
		  u=links3(mu).getReadable()),
	  A,
	  {
	    su3_copy(sm_conf[A][mu],u[A]);
	  });
  }
}
