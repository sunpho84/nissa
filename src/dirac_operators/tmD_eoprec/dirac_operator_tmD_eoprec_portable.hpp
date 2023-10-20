#ifndef _DIRAC_OPERATOR_TMD_EOPREC_PORTABLE_HPP
#define _DIRAC_OPERATOR_TMD_EOPREC_PORTABLE_HPP

#ifdef HAVE_CONFIG_H
# include "config.hpp"
#endif

#include "base/field.hpp"
#include "new_types/su3_op.hpp"

namespace nissa
{
  void tmDkern_eoprec_eos(OddField<spincolor>& out,
			  EvnField<spincolor>& temp,
			  const EoField<quad_su3>& conf,
			  const double& kappa,
			  const double& mu,
			  const OddField<spincolor>& in);
  
  void tmDkern_eoprec_square_eos(OddField<spincolor>& out,
				 OddField<spincolor>& temp1,
				 EvnField<spincolor>& temp2,
				 const EoField<quad_su3>& conf,
				 const double& kappa,
				 const double& mu,
				 const OddField<spincolor>& in);
  
  void tmDkern_eoprec_eos_put_together_and_include_gamma5(OddField<spincolor>& out,
							  const OddField<spincolor>& temp);
  
  /// Apply even-odd or odd-even part of tmD, multiplied by -2
  template <typename O,
	    typename I>
  void tmn2Deo_or_tmn2Doe_eos(O& out,
			      const EoField<quad_su3>& conf,
			      const I& in)
  {
    constexpr int xPar=O::fieldCoverage;
    static_assert(xPar!=I::fieldCoverage,"calling with messed up parities");
    
    conf.updateHalo();
    in.updateHalo();
    
    PAR(0,locVolh,
	CAPTURE(TO_WRITE(out),
		TO_READ(in),
		TO_READ(conf)),
	X,
	{
	  spincolor_put_to_zero(out[X]);
	  
	  for(int mu=0;mu<NDIM;mu++)
	    {
	      color0 temp_c0,temp_c1,temp_c2,temp_c3;
	      
	      //Forward
	      const int Xup=O::locNeighup(X,mu);
	      switch(mu)
		{
		case 0:
		  color_summ(temp_c0,in[Xup][0],in[Xup][2]);
		  color_summ(temp_c1,in[Xup][1],in[Xup][3]);
		  break;
		case 1:
		  color_isumm(temp_c0,in[Xup][0],in[Xup][3]);
		  color_isumm(temp_c1,in[Xup][1],in[Xup][2]);
		  break;
		case 2:
		  color_summ(temp_c0,in[Xup][0],in[Xup][3]);
		  color_subt(temp_c1,in[Xup][1],in[Xup][2]);
		  break;
		case 3:
		  color_isumm(temp_c0,in[Xup][0],in[Xup][2]);
		  color_isubt(temp_c1,in[Xup][1],in[Xup][3]);
		  break;
		}
	      
	      unsafe_su3_prod_color(temp_c2,conf[xPar][X][mu],temp_c0);
	      unsafe_su3_prod_color(temp_c3,conf[xPar][X][mu],temp_c1);
	      
	      color_summassign(out[X][0],temp_c2);
	      color_summassign(out[X][1],temp_c3);
	      
	      switch(mu)
		{
		case 0:
		  color_summassign(out[X][2],temp_c2);
		  color_summassign(out[X][3],temp_c3);
		  break;
		case 1:
		  color_isubtassign(out[X][2],temp_c3);
		  color_isubtassign(out[X][3],temp_c2);
		  break;
		case 2:
		  color_subtassign(out[X][2],temp_c3);
		  color_summassign(out[X][3],temp_c2);
		  break;
		case 3:
		  color_isubtassign(out[X][2],temp_c2);
		  color_isummassign(out[X][3],temp_c3);
		  break;
		}
	      
	      //Backward
	      const int Xdw=O::locNeighdw(X,mu);
	      switch(mu)
		{
		case 0:
		  color_subt(temp_c0,in[Xdw][0],in[Xdw][2]);
		  color_subt(temp_c1,in[Xdw][1],in[Xdw][3]);
		  break;
		case 1:
		  color_isubt(temp_c0,in[Xdw][0],in[Xdw][3]);
		  color_isubt(temp_c1,in[Xdw][1],in[Xdw][2]);
		  break;
		case 2:
		  color_subt(temp_c0,in[Xdw][0],in[Xdw][3]);
		  color_summ(temp_c1,in[Xdw][1],in[Xdw][2]);
		  break;
		case 3:
		  color_isubt(temp_c0,in[Xdw][0],in[Xdw][2]);
		  color_isumm(temp_c1,in[Xdw][1],in[Xdw][3]);
		  break;
		}
	      
	      unsafe_su3_dag_prod_color(temp_c2,conf[!xPar][Xdw][mu],temp_c0);
	      unsafe_su3_dag_prod_color(temp_c3,conf[!xPar][Xdw][mu],temp_c1);
	      
	      color_summassign(out[X][0],temp_c2);
	      color_summassign(out[X][1],temp_c3);
	      
	      switch(mu)
		{
		case 0:
		  color_subtassign(out[X][2],temp_c2);
		  color_subtassign(out[X][3],temp_c3);
		  break;
		case 1:
		  color_isummassign(out[X][2],temp_c3);
		  color_isummassign(out[X][3],temp_c2);
		  break;
		case 2:
		  color_summassign(out[X][2],temp_c3);
		  color_subtassign(out[X][3],temp_c2);
		  break;
		case 3:
		  color_isummassign(out[X][2],temp_c2);
		  color_isubtassign(out[X][3],temp_c3);
		  break;
		}
	    }
	});
  }
  
  //implement ee or oo part of Dirac operator, equation(3)
  template <typename O,
	    typename I>
  void tmDee_or_oo_eos(O&& out,
			   const double& kappa,
			   const double& mu,
			   const I& in)
  {
    PAR(0,locVolh,
	CAPTURE(kappa,mu,
		TO_WRITE(out),
		TO_READ(in)),
	X,
	{
	  for(int ic=0;ic<NCOL;ic++)
	    {
	      const complex z={1/(2*kappa),mu};
	      
	      for(int id=0;id<NDIRAC/2;id++)
		unsafe_complex_prod(out[X][id][ic],in[X][id][ic],z);
	      for(int id=NDIRAC/2;id<4;id++)
		unsafe_complex_conj2_prod(out[X][id][ic],in[X][id][ic],z);
	    }
	});
  }
  
  /// Inverse
  template <typename O,
	    typename I>
  void inv_tmDee_or_oo_eos(O&& out,
			   const double& kappa,
			   const double& mu,
			   const I& in)
  {
    PAR(0,locVolh,
	CAPTURE(kappa,mu,
		TO_WRITE(out),
		TO_READ(in)),
	X,
	{
	  for(int ic=0;ic<NCOL;ic++)
	    {
	      const double a=1/(2*kappa),b=mu,nrm=a*a+b*b;
	      const complex z={+a/nrm,-b/nrm};
	      
	      for(int id=0;id<NDIRAC/2;id++)
		unsafe_complex_prod(out[X][id][ic],in[X][id][ic],z);
	      for(int id=NDIRAC/2;id<4;id++)
		unsafe_complex_conj2_prod(out[X][id][ic],in[X][id][ic],z);
	    }
	});
  }
}

#endif
