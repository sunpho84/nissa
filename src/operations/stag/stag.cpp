#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "base/vectors.hpp"
#include "communicate/borders.hpp"
#include "geometry/geometry_eo.hpp"
#include "hmc/backfield.hpp"
#include "hmc/theory_pars.hpp"
#include "inverters/staggered/cg_invert_stD.hpp"
#include "linalgs/linalgs.hpp"
#include "new_types/su3.hpp"

#include "stag.hpp"

#ifdef USE_THREADS
 #include "routines/thread.hpp"
#endif

namespace nissa
{
  namespace stag
  {
    //multiply by M^-1
    THREADABLE_FUNCTION_6ARG(mult_Minv, color**,prop, quad_su3**,conf, quad_u1**,u1b, double,m, double,residue, color**,source)
    {
      add_backfield_to_conf(conf,u1b);
      inv_stD_cg(prop,conf,m,100000,residue,source);
      rem_backfield_from_conf(conf,u1b);
    }
    THREADABLE_FUNCTION_END
    void mult_Minv(color **prop,quad_su3 **conf,theory_pars_t *pars,int iflav,double residue,color **source)
    {mult_Minv(prop,conf,pars->backfield[iflav],pars->quarks[iflav].mass,residue,source);}
    
    //compute the matrix element of the derivative of the dirac operator between two vectors
    //forward and backward derivative are stored separately, for a reason
    void compute_fw_bw_der_mel(complex *res_fw_bw,color **left,quad_su3 **conf,int mu,color **right,complex *point_result)
    {
      GET_THREAD_ID();
      
      color **right_fw_bw[2]={right,left};
      
      for(int fw_bw=0;fw_bw<2;fw_bw++)
	{
	  vector_reset(point_result);
	  for(int par=0;par<2;par++)
	    NISSA_PARALLEL_LOOP(ieo,0,loc_volh)
	      {
		color v;
		unsafe_su3_prod_color(v,conf[par][ieo][mu],right_fw_bw[fw_bw][!par][loceo_neighup[par][ieo][mu]]);
		complex t;
		if(fw_bw==0) color_scalar_prod(t,right_fw_bw[!fw_bw][par][ieo],v);
		else         color_scalar_prod(t,v,right_fw_bw[!fw_bw][par][ieo]);
		complex_summassign(point_result[loclx_of_loceo[par][ieo]],t);
	      }
	  THREAD_BARRIER();
	  complex_vector_glb_collapse(res_fw_bw[fw_bw],point_result,loc_vol);
	}
    }
    
    //fill a source
    void fill_source(color **src)
    {generate_fully_undiluted_eo_source(src,RND_GAUSS,-1);}
    
    //take the trace between A^dag and B
    THREADABLE_FUNCTION_4ARG(summ_the_trace, double*,out, complex*,point_result, color**, A, color**, B)
    {
      GET_THREAD_ID();
      
      //compute results for single points
      vector_reset(point_result);
      for(int par=0;par<2;par++)
	NISSA_PARALLEL_LOOP(ieo,0,loc_volh)
	  for(int ic=0;ic<3;ic++)
	    complex_summ_the_conj1_prod(point_result[loclx_of_loceo[par][ieo]],A[par][ieo][ic],B[par][ieo][ic]);
      THREAD_BARRIER();
      
      //final reduction
      complex temp;
      complex_vector_glb_collapse(temp,point_result,loc_vol);
      if(IS_MASTER_THREAD) complex_summassign(out,temp);
    }
    THREADABLE_FUNCTION_END
    
    //multiply by the derivative of M w.r.t mu
    THREADABLE_FUNCTION_6ARG(mult_dMdmu, color**,out, theory_pars_t*,theory_pars, quad_su3**,conf, int,iflav, int,ord, color**,in)
    {
      GET_THREAD_ID();
      
      if(ord==0) crash("makes no sense to call with order zero");
      
      add_backfield_to_conf(conf,theory_pars->backfield[iflav]);
      communicate_ev_and_od_quad_su3_borders(conf);
      communicate_ev_and_od_color_borders(in);
      
      for(int par=0;par<2;par++)
	{
	  NISSA_PARALLEL_LOOP(ieo,0,loc_volh)
	    {
	      color temp;
	      unsafe_su3_prod_color(temp,conf[par][ieo][0],in[!par][loceo_neighup[par][ieo][0]]);
	      int idw=loceo_neighdw[par][ieo][0];
	      if(ord%2==0) su3_dag_subt_the_prod_color(temp,conf[!par][idw][0],in[!par][idw]);
	      else         su3_dag_summ_the_prod_color(temp,conf[!par][idw][0],in[!par][idw]);
	      color_prod_double(out[par][ieo],temp,0.5);
	    }
	  set_borders_invalid(out[par]);
	}
      
      rem_backfield_from_conf(conf,theory_pars->backfield[iflav]);
    }
    THREADABLE_FUNCTION_END
  
    void insert_external_source_handle(complex out,spin1field **aux,int par,int ieo,int mu,void *pars)
    {if(aux) complex_copy(out,aux[par][ieo][mu]);else complex_put_to_real(out,1);}
    //insert an external current
    void insert_vector_vertex(color **out,quad_su3 **conf,theory_pars_t *theory_pars,int iflav,spin1field **curr,color **in,complex fact_fw,complex fact_bw,void(*get_curr)(complex out,spin1field **curr,int par,int ieo,int mu,void *pars),int t,void *pars)
    {
      GET_THREAD_ID();
      
      add_backfield_to_conf(conf,theory_pars->backfield[iflav]);
      communicate_ev_and_od_quad_su3_borders(conf);
      communicate_ev_and_od_color_borders(in);
      
      for(int par=0;par<2;par++)
	{
	  NISSA_PARALLEL_LOOP(ieo,0,loc_volh)
	    {
	      color_put_to_zero(out[par][ieo]);
	      for(int mu=0;mu<NDIM;mu++)
		{
		  color temp;
		  
		  int iup=loceo_neighup[par][ieo][mu];
		  complex cf;
		  get_curr(cf,curr,par,ieo,mu,pars);
		  complex_prodassign(cf,fact_fw);
		  unsafe_su3_prod_color(temp,conf[par][ieo][mu],in[!par][iup]);
		  color_summ_the_prod_complex(out[par][ieo],temp,cf);
		  
		  int idw=loceo_neighdw[par][ieo][mu];
		  complex cb;
		  get_curr(cb,curr,!par,idw,mu,pars);
		  complex_prodassign(cb,fact_bw);
		  unsafe_su3_dag_prod_color(temp,conf[!par][idw][mu],in[!par][idw]);
		  color_summ_the_prod_complex(out[par][ieo],temp,cb);
		}
	    }
	  set_borders_invalid(out[par]);
	}
      
      rem_backfield_from_conf(conf,theory_pars->backfield[iflav]);
    }
  }
  
  std::string base_fermionic_meas_t::get_str(bool full)
  {
    std::ostringstream os;
    
    if(each!=def_each()||full) os<<" Each\t\t=\t"<<each<<"\n";
    if(after!=def_after()||full) os<<" After\t\t=\t"<<after<<"\n";
    if(path!=def_path()||full) os<<" Path\t\t=\t\""<<path.c_str()<<"\"\n";
    if(residue!=def_residue()||full) os<<" Residue\t=\t"<<residue<<"\n";
    if(ncopies!=def_ncopies()||full) os<<" NCopies\t=\t"<<ncopies<<"\n";
    if(itheory!=def_itheory()||full) os<<" ITheory\t=\t"<<itheory<<"\n";
    if(nhits!=def_nhits()||full) os<<" NHits\t\t=\t"<<nhits<<"\n";
    
    return os.str();
  }
}
