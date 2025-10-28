#ifndef _DIRAC_OPERATOR_STD_HPP
#define _DIRAC_OPERATOR_STD_HPP

#ifdef HAVE_CONFIG_H
# include "config.hpp"
#endif

#include "base/field.hpp"
#include "geometry/geometry_eo.hpp"
#include "new_types/su3_op.hpp"

namespace nissa
{
  /// Return the odd part of the application of D to a vector
  template <typename F>
  void odd_apply_stD(OddField<Color<F>>& out,
		     const EoField<QuadSu3<F>>& conf,
		     const double m,
		     const EoField<Color<F>>& in,
		     const double& sign=1)
  {
    apply_st2Doe(out,conf,in.evenPart);
    
    FOR_EACH_SITE_DEG_OF_FIELD(out,
			       CAPTURE(m,sign,
				       TO_WRITE(out),
				       TO_READ(in)),site,iDeg,
			       {
				 auto& d=out(site,iDeg);
				 d=in[ODD](site,iDeg)*m+0.5*sign*d;
			       });
  }
  
  template <typename F>
  void odd_apply_stD_dag(OddField<Color<F>>& out,
			 const EoField<QuadSu3<F>>& conf,
			 const double m,
			 const EoField<Color<F>>& in)
  {
    odd_apply_stD(out,conf,m,in,-1);
  }
  
  // void apply_Adams(eo_ptr<color> out,eo_ptr<quad_su3> conf,eo_ptr<quad_u1> u1b,double m,double m_Adams,eo_ptr<color> temp,eo_ptr<color> in);
  // void apply_AdamsII(eo_ptr<color> out,eo_ptr<quad_su3> conf,eo_ptr<quad_u1> u1b,double m,double m_Adams,eo_ptr<color> temp,eo_ptr<color> in);
  
  template <typename F>
  void apply_st2Doe(OddField<Color<F>>& out,
		    const EoField<QuadSu3<F>>& conf,
		    const EvnField<Color<F>>& in)
  {
    conf.updateHalo();
    in.updateHalo();
    
    PAR(0,locVolh,
	CAPTURE(TO_READ(conf),
		TO_READ(in),
		TO_WRITE(out)),io,
	{
	  //neighbours search
	  const int evup0=loceo_neighup[ODD][io][0];
	  const int evdw0=loceo_neighdw[ODD][io][0];
	  
	  //derivative in the time direction - without self-summ
	  unsafe_su3_prod_color(      out[io],conf[ODD][io   ][0],in[evup0]);
	  su3_dag_subt_the_prod_color(out[io],conf[EVN][evdw0][0],in[evdw0]);
	  
	  //derivatives in the spatial direction - with self summ
	  UNROLL_FOR(mu,1,NDIM)
	    {
	      const int evup=loceo_neighup[ODD][io][mu];
	      const int evdw=loceo_neighdw[ODD][io][mu];
	      
	      su3_summ_the_prod_color(    out[io],conf[ODD][io  ][mu],in[evup]);
	      su3_dag_subt_the_prod_color(out[io],conf[EVN][evdw][mu],in[evdw]);
	    }
	});
  }
  
  /// Put the 0.5 factor
  template <typename F>
  void apply_stDoe(OddField<Color<F>>& out,
		   const EoField<QuadSu3<F>>& conf,
		   const EvnField<Color<F>>& in)
  {
    apply_st2Doe(out,conf,in);
    
    PAR(0,locVolh,CAPTURE(TO_WRITE(out)),io,
	{
	  color_prodassign_double(out[io],0.5);
	});
  }
  
  template <typename F>
  void apply_stDeo_half(EvnField<Color<F>>& out,
			const EoField<QuadSu3<F>>& conf,
			const OddField<Color<F>>& in)
  {
    conf.updateHalo();
    in.updateHalo();
    
    PAR(0,locVolh,
	CAPTURE(TO_READ(conf),
		TO_WRITE(out),
		TO_READ(in)),ie,
	{
	  //neighbours search
	  const int odup0=loceo_neighup[EVN][ie][0];
	  const int oddw0=loceo_neighdw[EVN][ie][0];
	  
	  //derivative in the time direction - without self-summ
	  unsafe_su3_prod_color(      out[ie],conf[EVN][ie   ][0],in[odup0]);
	  su3_dag_subt_the_prod_color(out[ie],conf[ODD][oddw0][0],in[oddw0]);
	  
	  //derivatives in the spatial direction - with self summ
	  UNROLL_FOR(mu,1,NDIM)
	    {
	      const int odup=loceo_neighup[EVN][ie][mu];
	      const int oddw=loceo_neighdw[EVN][ie][mu];
	      
	      su3_summ_the_prod_color(    out[ie],conf[EVN][ie  ][mu],in[odup]);
	      su3_dag_subt_the_prod_color(out[ie],conf[ODD][oddw][mu],in[oddw]);
	    }
	  
	  //Doe contains 1/2, we put an additional one
	  color_prod_double(out[ie],out[ie],0.25);
	});
  }
  
  template <typename F>
  void apply_stD2ee_m2(EvnField<Color<F>>& out,
		       const EoField<QuadSu3<F>>& conf,
		       OddField<Color<F>>& temp,
		       const double& mass2,
		       const EvnField<Color<F>>& in)
  {
    //check arguments
    if(out==in)   CRASH("out==in!");
    START_TIMING(portable_stD_app_time,nportable_stD_app);
    
    conf.updateHalo();
    in.updateHalo();
    
    PAR(0,locVolh,
	CAPTURE(TO_READ(in),
		TO_READ(conf),
		TO_WRITE(temp)),io,
	{
	  Color<F> o{};
	  
	  //derivatives in the spatial direction - with self summ
	  UNROLL_FOR_ALL_DIRS(mu)
	    {
	      const int evup=loceo_neighup[ODD][io][mu];
	      const int evdw=loceo_neighdw[ODD][io][mu];
	      
	      Color<F> iu,id;
	      color_copy(iu,in[evup]);
	      color_copy(id,in[evdw]);
	      
	      Su3<F> cu,cd;
	      su3_copy(cu,conf[ODD][io  ][mu]);
	      su3_copy(cd,conf[EVN][evdw][mu]);
	      
	      su3_summ_the_prod_color(o,cu,iu);
	      su3_dag_subt_the_prod_color(o,cd,id);
	    }
	  
	  color_copy(temp[io],o);
	});
    
    temp.updateHalo();
    
    //we still apply Deo, but then we put a - because we should apply Doe^+=-Deo
    PAR(0,locVolh,
	CAPTURE(TO_READ(temp),
		TO_READ(conf),
		TO_WRITE(out)),ie,
	{
	  Color<F>o{};
	  
	  UNROLL_FOR_ALL_DIRS(mu)
	    {
	      const int odup=loceo_neighup[EVN][ie][mu];
	      const int oddw=loceo_neighdw[EVN][ie][mu];
	      
	      Color<F>tu,td;
	      color_copy(tu,temp[odup]);
	      color_copy(td,temp[oddw]);
	      
	      Su3<F> cu,cd;
	      su3_copy(cu,conf[EVN][ie  ][mu]);
	      su3_copy(cd,conf[ODD][oddw][mu]);
	      
	      su3_summ_the_prod_color(o,cu,tu);
	      su3_dag_subt_the_prod_color(o,cd,td);
	    }
	  
	  color_copy(out[ie],o);
	});
    
    if(mass2!=0)
      {
	PAR(0,locVolh,
	    CAPTURE(mass2,
		    TO_WRITE(out),
		    TO_READ(in)),ie,
	    {
	      UNROLL_FOR_ALL_COLS(ic)
		UNROLL_FOR(ri,0,2)
		  out[ie][ic][ri]=mass2*in[ie][ic][ri]-out[ie][ic][ri]*0.25;
	    });
      }
    else
      {
	PAR(0,locVolh,
	    CAPTURE(TO_WRITE(out)),ie,
	    {
	      UNROLL_FOR_ALL_COLS(ic)
		UNROLL_FOR(ri,0,2)
		  out[ie][ic][ri]*=-0.25;
	    });
      }
    
    STOP_TIMING(portable_stD_app_time);
  }
  
  /// Functor needed to find the call to apply_stD2ee_m2 with specific conf and temp
  struct ApplyStD2eeM2Functor
  {
    /// Conf to be used
    const EoField<quad_su3>& eo_conf;
    
    /// Temporary vector
    OddField<color>& temp;
    
    /// Callable
    void operator()(EvnField<color>& out,
		    const double& mass2,
		    const EvnField<color>& in)
    {
      apply_stD2ee_m2(out,eo_conf,temp,mass2,in);
    }
    
    /// Construct taking reference
    ApplyStD2eeM2Functor(const EoField<quad_su3>& eo_conf,
			 OddField<color>& temp) :
      eo_conf(eo_conf),
      temp(temp)
    {
    }
  };
  
  /// Return the even part of the application of D to a vector
  template <typename F>
  void evn_apply_stD(EvnField<Color<F>>& out,
		     const EoField<QuadSu3<F>>& conf,
		     const double m,
		     const EoField<Color<F>>& in,
		     const double sign=1)
  {
    apply_stDeo_half(out,conf,in.oddPart);
    
    FOR_EACH_SITE_DEG_OF_FIELD(out,
			       CAPTURE(m,sign,
				       TO_WRITE(out),
				       TO_READ(in)),site,iDeg,
			       {
				 auto& d= out(site,iDeg);
				 d=in[EVN](site,iDeg)*m+2*sign*d;
			       });
  }
  
  template <typename F>
  void evn_apply_stD_dag(EvnField<Color<F>>& out,
			 const EoField<QuadSu3<F>>& conf,
			 const double m,
			 const EoField<Color<F>>& in)
  {
    evn_apply_stD(out,conf,m,in,-1);
  }
  
  /// Return the result of the application of D to a vector
  template <typename F>
  void apply_stD(EoField<Color<F>>& out,
		 const EoField<QuadSu3<F>>& conf,
		 const double& m,
		 const EoField<Color<F>>& in)
  {
    evn_apply_stD(out.evenPart,conf,m,in);
    odd_apply_stD(out.oddPart,conf,m,in);
  }
}

#endif
