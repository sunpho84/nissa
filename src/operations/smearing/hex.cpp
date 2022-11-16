#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

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
  template <typename LIn,
	    typename...Tail>
  void sme(su3* lOut,const double& alpha,const LIn& V,su3* U,const Dir& mu,const Tail&...a)
  {
    vector_reset(lOut);
    
    for(Dir sigma=0;sigma<NDIM;sigma++)
      {
	bool toDo=true;
	for(const auto& i : {mu,a...})
	  toDo&=(i!=sigma);
	
	const auto& A=V(sigma,mu,a...);
	const auto& B=V(mu,a...,sigma);
	
	if(toDo)
	  NISSA_PARALLEL_LOOP(iVol,0,locVol)
	    {
	      su3 temp;
	      
	      const LocLxSite muUp=loclxNeighup(iVol,mu);
	      
	      unsafe_su3_prod_su3(temp,A[iVol],B[loclxNeighup(iVol,sigma).nastyConvert()]);
	      su3_summ_the_prod_su3_dag(lOut[iVol.nastyConvert()],temp,A[muUp]);
	      
	      const LocLxSite sigmaDw=loclxNeighdw(iVol,sigma);
	      const LocLxSite sigmaDwMuUp=loclxNeighdw(muUp,sigma);
	      
	      unsafe_su3_dag_prod_su3(temp,A[sigmaDw],B[sigmaDw]);
	      su3_summ_the_prod_su3(lOut[iVol.nastyConvert()],temp,A[sigmaDwMuUp]);
	    }
	NISSA_PARALLEL_LOOP_END;
      }
    
    NISSA_PARALLEL_LOOP(iVol,0,locVol)
      {
	su3 temp1,temp2;
	
	unsafe_su3_prod_su3_dag(temp1,lOut[iVol.nastyConvert()],U[iVol.nastyConvert()]);
	unsafe_su3_traceless_anti_hermitian_part(temp2,temp1);
	
	su3 Q;
	unsafe_su3_traceless_anti_hermitian_part(Q,temp2);
	su3_prodassign_idouble(Q,-alpha);
	
	safe_hermitian_exact_i_exponentiate(Q,Q);
	unsafe_su3_prod_su3(lOut[iVol.nastyConvert()],Q,U[iVol.nastyConvert()]);
      }
    NISSA_PARALLEL_LOOP_END;
  }
  
  /// Hex smear the conf
  void hex_smear_conf(quad_su3* sm_conf,quad_su3* conf,const double& alpha1,const double& alpha2,const double& alpha3)
  {
    crash("must be done inside a struct");
    
    // const int nDecLevels=4;
    // const int nDecsPerLevel[]=
    //   {NDIM,
    //    NDIM*(NDIM-1),
    //    NDIM*(NDIM-1),
    //    NDIM};
    
    // const int nDecTot=
    // 		nDecsPerLevel[0]+nDecsPerLevel[1]+nDecsPerLevel[2]+nDecsPerLevel[3];
    // su3* linksAllDecs[nDecTot];
    
    // for(int iDec=0;iDec<nDecTot;iDec++)
    //   linksAllDecs[iDec]=nissa_malloc("Dec",locVol+bord_vol+edge_vol,su3);
    
    // su3** linksDecsPerLev[nDecLevels];
    // linksDecsPerLev[0]=linksAllDecs;
    // for(int iDecLev=1;iDecLev<nDecLevels;iDecLev++)
    //   linksDecsPerLev[iDecLev]=linksDecsPerLev[iDecLev-1]+nDecsPerLevel[iDecLev-1];
    
    // auto communicateDecsForLev=
    //   [linksDecsPerLev,nDecsPerLevel](const int& iLev)
    //   {
    // 	for(int iDec=0;iDec<nDecsPerLevel[iLev];iDec++)
    // 	  {
    // 	    su3* o=linksDecsPerLev[iLev][iDec];
    // 	    set_borders_invalid(o);
    // 	    communicate_lx_su3_edges(o);
    // 	  }
    //   };
    
    // /////////////////////////////////////////////////////////////////
    
    // /// Level0
    
    // const auto links0=
    //   [links0List=linksDecsPerLev[0]](const int& mu,auto...) //swallow all other dirs
    //   {
    // 	return
    // 	  links0List[mu];
    //   };
    
    // NISSA_PARALLEL_LOOP(A,0,locVol)
    //   for(int mu=0;mu<NDIM;mu++)
    // 	su3_copy(links0(mu)[A],conf[A][mu]);
    // NISSA_PARALLEL_LOOP_END;
    
    // communicateDecsForLev(0);
    
    // /////////////////////////////////////////////////////////////////
    
    // /// Level1
    
    // const auto links1=
    //   [links1List=linksDecsPerLev[1]](const int& mu,const int& nu,const int& rho)
    //   {
    // 	const int e=15^((1<<mu)+(1<<nu)+(1<<rho));
    // 	const int sigma=(e==2)+2*(e==4)+3*(e==8);
    // 	const int iSigma=sigma-(sigma>mu); /// Subtract 1 if sigma>mu since we don't store mu==sigma
	
    // 	return
    // 	  links1List[iSigma+(NDIM-1)*mu];
    //   };
    
    // for(int mu=0;mu<NDIM;mu++)
    //   for(int iNu=0;iNu<NDIM-1;iNu++)
    // 	for(int iRho=0;iRho<NDIM-2;iRho++)
    // 	  {
    // 	    const int nu=perp_dir[mu][iNu];
    // 	    const int rho=perp2_dir[mu][iNu][iRho];
	    
    // 	    if(rho>nu)
    // 	      sme(links1(mu,nu,rho),alpha3/2,links0,links0(mu),mu,nu,rho);
    // 	  }
    
    // communicateDecsForLev(1);
    
    // /////////////////////////////////////////////////////////////////
    
    // /// Level2
    
    // const auto links2=
    //   [links2List=linksDecsPerLev[2]](const int& mu,const int& nu)
    //   {
    // 	const int iNu=
    // 	  nu-(nu>mu); /// Subtract 1 if sigma>mu since we don't store mu==sigma
	
    // 	return
    // 	  links2List[iNu+(NDIM-1)*mu];
    //   };
    
    // for(int mu=0;mu<NDIM;mu++)
    //   for(int iNu=0;iNu<NDIM-1;iNu++)
    // 	{
    // 	  const int nu=perp_dir[mu][iNu];
    // 	  sme(links2(mu,nu),alpha2/4,links1,links0(mu),mu,nu);
    // 	}
    
    // communicateDecsForLev(2);
    
    // /////////////////////////////////////////////////////////////////
    
    // /// Level3
    
    // const auto links3=
    //   [links3List=linksDecsPerLev[3]](const int& mu)
    //   {
    // 	return
    // 	  links3List[mu];
    //   };
    
    // for(int mu=0;mu<NDIM;mu++)
    //   sme(links3(mu),alpha1/6,links2,links0(mu),mu);
    
    // /////////////////////////////////////////////////////////////////
    
    // for(int mu=0;mu<NDIM;mu++)
    //   NISSA_PARALLEL_LOOP(A,0,locVol)
    // 	su3_copy(sm_conf[A][mu],links3(mu)[A]);
    // NISSA_PARALLEL_LOOP_END;
    
    // set_borders_invalid(sm_conf);
    
    // for(int iDec=0;iDec<nDecTot;iDec++)
    //   nissa_free(linksAllDecs[iDec]);
  }
}
