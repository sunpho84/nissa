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
  //apply even-odd or odd-even part of tmD, multiplied by -2
  void tmn2Deo_or_tmn2Doe_eos(spin *out,int eooe,spin *in,momentum_t bc)
  {
    GET_THREAD_ID();
    
    if(eooe==0) communicate_od_spin_borders(in);
    else        communicate_ev_spin_borders(in);
    
    int Xup,Xdw;
    
    complex phases[4];
    for(int mu=0;mu<4;mu++)
      {    
	phases[mu][RE]=cos(M_PI*bc[mu]);
	phases[mu][IM]=sin(M_PI*bc[mu]);
      }
    
    NISSA_PARALLEL_LOOP(X,0,loc_volh)
      {
	complex temp_c0,temp_c1;
	
	//Forward 0
	Xup=loceo_neighup[eooe][X][0];
	complex_summ(temp_c0,in[Xup][0],in[Xup][2]);
	complex_summ(temp_c1,in[Xup][1],in[Xup][3]);
	if(glb_coord_of_loclx[loclx_of_loceo[!eooe][Xup]][0]==0)
	  {
	    safe_complex_prod(temp_c0,temp_c0,phases[0]);
	    safe_complex_prod(temp_c1,temp_c1,phases[0]);
	  }
	complex_copy(out[X][0],temp_c0);
	complex_copy(out[X][1],temp_c1);
	complex_copy(out[X][2],out[X][0]);
	complex_copy(out[X][3],out[X][1]);
	
	//Backward 0
	Xdw=loceo_neighdw[eooe][X][0];
	complex_subt(temp_c0,in[Xdw][0],in[Xdw][2]);
	complex_subt(temp_c1,in[Xdw][1],in[Xdw][3]);
	if(glb_coord_of_loclx[loclx_of_loceo[eooe][X]][0]==0)
	  {
	    safe_complex_conj2_prod(temp_c0,temp_c0,phases[0]);
	    safe_complex_conj2_prod(temp_c1,temp_c1,phases[0]);
	  }
	complex_summassign(out[X][0],temp_c0);
	complex_summassign(out[X][1],temp_c1);
	complex_subtassign(out[X][2],temp_c0);
	complex_subtassign(out[X][3],temp_c1);
	
	//Forward 1
	Xup=loceo_neighup[eooe][X][1];
	complex_isumm(temp_c0,in[Xup][0],in[Xup][3]);
	complex_isumm(temp_c1,in[Xup][1],in[Xup][2]);
	if(glb_coord_of_loclx[loclx_of_loceo[!eooe][Xup]][1]==0)
	  {
	    safe_complex_prod(temp_c0,temp_c0,phases[1]);
	    safe_complex_prod(temp_c1,temp_c1,phases[1]);
	  }
	complex_summassign(out[X][0],temp_c0);
	complex_summassign(out[X][1],temp_c1);
	complex_isubtassign(out[X][2],temp_c1);
	complex_isubtassign(out[X][3],temp_c0);
	
	//Backward 1
	Xdw=loceo_neighdw[eooe][X][1];
	complex_isubt(temp_c0,in[Xdw][0],in[Xdw][3]);
	complex_isubt(temp_c1,in[Xdw][1],in[Xdw][2]);
	if(glb_coord_of_loclx[loclx_of_loceo[eooe][X]][1]==0)
	  {
	    safe_complex_conj2_prod(temp_c0,temp_c0,phases[1]);
	    safe_complex_conj2_prod(temp_c1,temp_c1,phases[1]);
	  }
	complex_summassign(out[X][0],temp_c0);
	complex_summassign(out[X][1],temp_c1);
	complex_isummassign(out[X][2],temp_c1);
	complex_isummassign(out[X][3],temp_c0);
	
	//Forward 2
	Xup=loceo_neighup[eooe][X][2];
	complex_summ(temp_c0,in[Xup][0],in[Xup][3]);
	complex_subt(temp_c1,in[Xup][1],in[Xup][2]);
	if(glb_coord_of_loclx[loclx_of_loceo[!eooe][Xup]][2]==0)
	  {
	    safe_complex_prod(temp_c0,temp_c0,phases[2]);
	    safe_complex_prod(temp_c1,temp_c1,phases[2]);
	  }
	complex_summassign(out[X][0],temp_c0);
	complex_summassign(out[X][1],temp_c1);
	complex_subtassign(out[X][2],temp_c1);
	complex_summassign(out[X][3],temp_c0);
	
	//Backward 2
	Xdw=loceo_neighdw[eooe][X][2];
	complex_subt(temp_c0,in[Xdw][0],in[Xdw][3]);
	complex_summ(temp_c1,in[Xdw][1],in[Xdw][2]);
	if(glb_coord_of_loclx[loclx_of_loceo[eooe][X]][2]==0)
	  {
	    safe_complex_conj2_prod(temp_c0,temp_c0,phases[2]);
	    safe_complex_conj2_prod(temp_c1,temp_c1,phases[2]);
	  }
	complex_summassign(out[X][0],temp_c0);
	complex_summassign(out[X][1],temp_c1);
	complex_summassign(out[X][2],temp_c1);
	complex_subtassign(out[X][3],temp_c0);
	
	//Forward 3
	Xup=loceo_neighup[eooe][X][3];
	complex_isumm(temp_c0,in[Xup][0],in[Xup][2]);
	complex_isubt(temp_c1,in[Xup][1],in[Xup][3]);
	if(glb_coord_of_loclx[loclx_of_loceo[!eooe][Xup]][3]==0)
	  {
	    safe_complex_prod(temp_c0,temp_c0,phases[3]);
	    safe_complex_prod(temp_c1,temp_c1,phases[3]);
	  }
	complex_summassign(out[X][0],temp_c0);
	complex_summassign(out[X][1],temp_c1);
	complex_isubtassign(out[X][2],temp_c0);
	complex_isummassign(out[X][3],temp_c1);
	
	//Backward 3
	Xdw=loceo_neighdw[eooe][X][3];
	complex_isubt(temp_c0,in[Xdw][0],in[Xdw][2]);
	complex_isumm(temp_c1,in[Xdw][1],in[Xdw][3]);
	if(glb_coord_of_loclx[loclx_of_loceo[eooe][X]][3]==0)
	  {
	    safe_complex_conj2_prod(temp_c0,temp_c0,phases[3]);
	    safe_complex_conj2_prod(temp_c1,temp_c1,phases[3]);
	  }
	complex_summassign(out[X][0],temp_c0);
	complex_summassign(out[X][1],temp_c1);
	complex_isummassign(out[X][2],temp_c0);
	complex_isubtassign(out[X][3],temp_c1);
      }
    NISSA_PARALLEL_LOOP_END;
    
    set_borders_invalid(out);
  }
  
  //wrappers
  void tmn2Doe_eos(spin *out,spin *in,momentum_t bc){tmn2Deo_or_tmn2Doe_eos(out,1,in,bc);}
  void tmn2Deo_eos(spin *out,spin *in,momentum_t bc){tmn2Deo_or_tmn2Doe_eos(out,0,in,bc);}
  
  //implement ee or oo part of Dirac operator, equation(3)
  void tmDee_or_oo_eos(spin *out,tm_quark_info qu,spin *in)
  {
    GET_THREAD_ID();
    
    if(in==out) crash("in==out!");
    complex z={1/(2*qu.kappa),qu.mass};
    
    NISSA_PARALLEL_LOOP(X,0,loc_volh)
      {
	for(int id=0;id<2;id++) unsafe_complex_prod(out[X][id],in[X][id],z);
	for(int id=2;id<4;id++) unsafe_complex_conj2_prod(out[X][id],in[X][id],z);
      }
    NISSA_PARALLEL_LOOP_END;
    
    set_borders_invalid(out);
  }
  
  //inverse
  void inv_tmDee_or_oo_eos(spin *out,tm_quark_info qu,spin *in)
  {
    GET_THREAD_ID();
    
    if(in==out) crash("in==out!");
    double a=1/(2*qu.kappa),b=qu.mass,nrm=a*a+b*b;
    complex z={+a/nrm,-b/nrm};
    
    NISSA_PARALLEL_LOOP(X,0,loc_volh)
      {
	for(int id=0;id<2;id++) unsafe_complex_prod(out[X][id],in[X][id],z);
	for(int id=2;id<4;id++) unsafe_complex_conj2_prod(out[X][id],in[X][id],z);
      }
    NISSA_PARALLEL_LOOP_END;
    
    set_borders_invalid(out);
  }
  
  //implement Koo defined in equation (7) 
  void tmDkern_eoprec_eos(spin *out,spin *temp,tm_quark_info qu,spin *in)
  {
    GET_THREAD_ID();
    
    tmn2Deo_eos(out,in,qu.bc);
    inv_tmDee_or_oo_eos(temp,qu,out);
    tmn2Doe_eos(out,temp,qu.bc);
    inv_tmDee_or_oo_eos(temp,qu,out);
    tmDee_or_oo_eos(temp,qu,in);
    
    NISSA_PARALLEL_LOOP(ivol,0,loc_volh)
      for(int id=0;id<2;id++)
	for(int ri=0;ri<2;ri++)
	  { //gamma5 is explicitely implemented
	    out[ivol][id  ][ri]=+temp[ivol][id  ][ri]-out[ivol][id  ][ri]*0.25;
	    out[ivol][id+2][ri]=-temp[ivol][id+2][ri]+out[ivol][id+2][ri]*0.25;
	  }
    NISSA_PARALLEL_LOOP_END;
    
    set_borders_invalid(out);
  }
  
  //square of Koo
  void tmDkern_eoprec_square_eos(spin *out,spin *temp1,spin *temp2,tm_quark_info qu,spin *in)
  {
    tm_quark_info mqu=qu;
    mqu.mass*=-1;
    
    tmDkern_eoprec_eos(temp1,temp2,mqu, in   );
    tmDkern_eoprec_eos(out,  temp2,qu,  temp1);
  }
}
