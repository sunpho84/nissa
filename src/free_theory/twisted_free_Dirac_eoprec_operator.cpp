#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <math.h>

#include "base/debug.hpp"
#include "base/vectors.hpp"
#include "communicate/borders.hpp"
#include "geometry/geometry_eo.hpp"
#include "new_types/complex.hpp"
#include "threads/threads.hpp"

#include "free_theory_types.hpp"

namespace nissa
{
  void tmn2Deo_or_tmn2Doe_eos(spin *out,const Parity& eooe,spin *in,const Momentum& bc)
  {
    if(eooe==0) communicate_od_spin_borders(in);
    else        communicate_ev_spin_borders(in);
    
    std::array<complex,4> phases;
    FOR_ALL_DIRS(mu)
      {
	phases[mu.nastyConvert()][RE]=cos(M_PI*bc(mu));
	phases[mu.nastyConvert()][IM]=sin(M_PI*bc(mu));
      }
    
    NISSA_PARALLEL_LOOP(X,0,locVolh)
      {
	LocEoSite Xup,Xdw;
	
	complex temp_c0,temp_c1;
	
	//Forward 0
	Xup=loceo_neighup(eooe,X,Dir(0));
	complex_summ(temp_c0,in[Xup.nastyConvert()][0],in[Xup.nastyConvert()][2]);
	complex_summ(temp_c1,in[Xup.nastyConvert()][1],in[Xup.nastyConvert()][3]);
	if(glbCoordOfLoclx(loclx_of_loceo(1-eooe,Xup),Dir(0))==0)
	  {
	    safe_complex_prod(temp_c0,temp_c0,phases[0]);
	    safe_complex_prod(temp_c1,temp_c1,phases[0]);
	  }
	complex_copy(out[X.nastyConvert()][0],temp_c0);
	complex_copy(out[X.nastyConvert()][1],temp_c1);
	complex_copy(out[X.nastyConvert()][2],out[X.nastyConvert()][0]);
	complex_copy(out[X.nastyConvert()][3],out[X.nastyConvert()][1]);
	
	//Backward 0
	Xdw=loceo_neighdw(eooe,X,Dir(0));
	complex_subt(temp_c0,in[Xdw.nastyConvert()][0],in[Xdw.nastyConvert()][2]);
	complex_subt(temp_c1,in[Xdw.nastyConvert()][1],in[Xdw.nastyConvert()][3]);
	if(glbCoordOfLoclx(loclx_of_loceo(eooe,X),Dir(0))==0)
	  {
	    safe_complex_conj2_prod(temp_c0,temp_c0,phases[0]);
	    safe_complex_conj2_prod(temp_c1,temp_c1,phases[0]);
	  }
	complex_summassign(out[X.nastyConvert()][0],temp_c0);
	complex_summassign(out[X.nastyConvert()][1],temp_c1);
	complex_subtassign(out[X.nastyConvert()][2],temp_c0);
	complex_subtassign(out[X.nastyConvert()][3],temp_c1);
	
	//Forward 1
	Xup=loceo_neighup(eooe,X,xDir);
	complex_isumm(temp_c0,in[Xup.nastyConvert()][0],in[Xup.nastyConvert()][3]);
	complex_isumm(temp_c1,in[Xup.nastyConvert()][1],in[Xup.nastyConvert()][2]);
	if(glbCoordOfLoclx(loclx_of_loceo(1-eooe,Xup),xDir)==0)
	  {
	    safe_complex_prod(temp_c0,temp_c0,phases[1]);
	    safe_complex_prod(temp_c1,temp_c1,phases[1]);
	  }
	complex_summassign(out[X.nastyConvert()][0],temp_c0);
	complex_summassign(out[X.nastyConvert()][1],temp_c1);
	complex_isubtassign(out[X.nastyConvert()][2],temp_c1);
	complex_isubtassign(out[X.nastyConvert()][3],temp_c0);
	
	//Backward 1
	Xdw=loceo_neighdw(eooe,X,xDir);
	complex_isubt(temp_c0,in[Xdw.nastyConvert()][0],in[Xdw.nastyConvert()][3]);
	complex_isubt(temp_c1,in[Xdw.nastyConvert()][1],in[Xdw.nastyConvert()][2]);
	if(glbCoordOfLoclx(loclx_of_loceo(eooe,X),xDir)==0)
	  {
	    safe_complex_conj2_prod(temp_c0,temp_c0,phases[1]);
	    safe_complex_conj2_prod(temp_c1,temp_c1,phases[1]);
	  }
	complex_summassign(out[X.nastyConvert()][0],temp_c0);
	complex_summassign(out[X.nastyConvert()][1],temp_c1);
	complex_isummassign(out[X.nastyConvert()][2],temp_c1);
	complex_isummassign(out[X.nastyConvert()][3],temp_c0);
	
	//Forward 2
	Xup=loceo_neighup(eooe,X,yDir);
	complex_summ(temp_c0,in[Xup.nastyConvert()][0],in[Xup.nastyConvert()][3]);
	complex_subt(temp_c1,in[Xup.nastyConvert()][1],in[Xup.nastyConvert()][2]);
	if(glbCoordOfLoclx(loclx_of_loceo(1-eooe,X),yDir)==0)
	  {
	    safe_complex_prod(temp_c0,temp_c0,phases[2]);
	    safe_complex_prod(temp_c1,temp_c1,phases[2]);
	  }
	complex_summassign(out[X.nastyConvert()][0],temp_c0);
	complex_summassign(out[X.nastyConvert()][1],temp_c1);
	complex_subtassign(out[X.nastyConvert()][2],temp_c1);
	complex_summassign(out[X.nastyConvert()][3],temp_c0);
	
	//Backward 2
	Xdw=loceo_neighdw(eooe,X,yDir);
	complex_subt(temp_c0,in[Xdw.nastyConvert()][0],in[Xdw.nastyConvert()][3]);
	complex_summ(temp_c1,in[Xdw.nastyConvert()][1],in[Xdw.nastyConvert()][2]);
	if(glbCoordOfLoclx(loclx_of_loceo(eooe,X),yDir)==0)
	  {
	    safe_complex_conj2_prod(temp_c0,temp_c0,phases[2]);
	    safe_complex_conj2_prod(temp_c1,temp_c1,phases[2]);
	  }
	complex_summassign(out[X.nastyConvert()][0],temp_c0);
	complex_summassign(out[X.nastyConvert()][1],temp_c1);
	complex_summassign(out[X.nastyConvert()][2],temp_c1);
	complex_subtassign(out[X.nastyConvert()][3],temp_c0);
	
	//Forward 3
	Xup=loceo_neighup(eooe,X,zDir);
	complex_isumm(temp_c0,in[Xup.nastyConvert()][0],in[Xup.nastyConvert()][2]);
	complex_isubt(temp_c1,in[Xup.nastyConvert()][1],in[Xup.nastyConvert()][3]);
	if(glbCoordOfLoclx(loclx_of_loceo(1-eooe,X),zDir)==0)
	  {
	    safe_complex_prod(temp_c0,temp_c0,phases[3]);
	    safe_complex_prod(temp_c1,temp_c1,phases[3]);
	  }
	complex_summassign(out[X.nastyConvert()][0],temp_c0);
	complex_summassign(out[X.nastyConvert()][1],temp_c1);
	complex_isubtassign(out[X.nastyConvert()][2],temp_c0);
	complex_isummassign(out[X.nastyConvert()][3],temp_c1);
	
	//Backward 3
	Xdw=loceo_neighdw(eooe,X,zDir);
	complex_isubt(temp_c0,in[Xdw.nastyConvert()][0],in[Xdw.nastyConvert()][2]);
	complex_isumm(temp_c1,in[Xdw.nastyConvert()][1],in[Xdw.nastyConvert()][3]);
	if(glbCoordOfLoclx(loclx_of_loceo(eooe,X),zDir)==0)
	  {
	    safe_complex_conj2_prod(temp_c0,temp_c0,phases[3]);
	    safe_complex_conj2_prod(temp_c1,temp_c1,phases[3]);
	  }
	complex_summassign(out[X.nastyConvert()][0],temp_c0);
	complex_summassign(out[X.nastyConvert()][1],temp_c1);
	complex_isummassign(out[X.nastyConvert()][2],temp_c0);
	complex_isubtassign(out[X.nastyConvert()][3],temp_c1);
      }
    NISSA_PARALLEL_LOOP_END;
    
    set_borders_invalid(out);
  }
  
  //wrappers
  void tmn2Doe_eos(spin *out,spin *in,const Momentum& bc)
  {
    tmn2Deo_or_tmn2Doe_eos(out,1,in,bc);
  }
  
  void tmn2Deo_eos(spin *out,spin *in,const Momentum& bc)
  {
    tmn2Deo_or_tmn2Doe_eos(out,0,in,bc);
  }
  
  //implement ee or oo part of Dirac operator, equation(3)
  void tmDee_or_oo_eos(spin *out,const tm_quark_info& qu,spin *in)
  {
    if(in==out) crash("in==out!");
    
    NISSA_PARALLEL_LOOP(X,0,locVolh)
      {
	const complex z={1/(2*qu.kappa),qu.mass};
	
	for(int id=0;id<2;id++) unsafe_complex_prod(out[X.nastyConvert()][id],in[X.nastyConvert()][id],z);
	for(int id=2;id<4;id++) unsafe_complex_conj2_prod(out[X.nastyConvert()][id],in[X.nastyConvert()][id],z);
      }
    NISSA_PARALLEL_LOOP_END;
    
    set_borders_invalid(out);
  }
  
  //inverse
  void inv_tmDee_or_oo_eos(spin *out,const tm_quark_info& qu,spin *in)
  {
    
    if(in==out) crash("in==out!");
    const double a=1/(2*qu.kappa),b=qu.mass,nrm=1/(a*a+b*b);
    
    NISSA_PARALLEL_LOOP(X,0,locVolh)
      {
	const complex z={+a*nrm,-b*nrm};
	for(int id=0;id<2;id++) unsafe_complex_prod(out[X.nastyConvert()][id],in[X.nastyConvert()][id],z);
	for(int id=2;id<4;id++) unsafe_complex_conj2_prod(out[X.nastyConvert()][id],in[X.nastyConvert()][id],z);
      }
    NISSA_PARALLEL_LOOP_END;
    
    set_borders_invalid(out);
  }
  
  //implement Koo defined in equation (7) 
  void tmDkern_eoprec_eos(spin *out,spin *temp,const tm_quark_info& qu,spin *in)
  {
    
    tmn2Deo_eos(out,in,qu.bc);
    inv_tmDee_or_oo_eos(temp,qu,out);
    tmn2Doe_eos(out,temp,qu.bc);
    inv_tmDee_or_oo_eos(temp,qu,out);
    tmDee_or_oo_eos(temp,qu,in);
    
    NISSA_PARALLEL_LOOP(ieo,0,locVolh)
      for(int id=0;id<2;id++)
	for(int ri=0;ri<2;ri++)
	  {
	    //gamma5 is explicitely implemented
	    out[ieo.nastyConvert()][id  ][ri]=+temp[ieo.nastyConvert()][id  ][ri]-out[ieo.nastyConvert()][id  ][ri]*0.25;
	    out[ieo.nastyConvert()][id+2][ri]=-temp[ieo.nastyConvert()][id+2][ri]+out[ieo.nastyConvert()][id+2][ri]*0.25;
	  }
    NISSA_PARALLEL_LOOP_END;
    
    set_borders_invalid(out);
  }
  
  //square of Koo
  void tmDkern_eoprec_square_eos(spin *out,spin *temp1,spin *temp2,const tm_quark_info& qu,spin *in)
  {
    tm_quark_info mqu;
    mqu.nastyCopy(qu);
    mqu.mass*=-1;
    
    tmDkern_eoprec_eos(temp1,temp2,mqu, in   );
    tmDkern_eoprec_eos(out,  temp2,qu,  temp1);
  }
}
