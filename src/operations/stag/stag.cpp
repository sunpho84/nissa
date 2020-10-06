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
    //summ a single shift
    void summ_covariant_shift(color **out,quad_su3 **conf,int mu,color **in,shift_orie_t side)
    {
      GET_THREAD_ID();
      
      if(in==out) crash("in==out");
      
      communicate_ev_and_od_color_borders(in);
      communicate_ev_and_od_quad_su3_borders(conf);
      
      for(int eo=0;eo<2;eo++)
	{
	  NISSA_PARALLEL_LOOP(ieo,0,loc_volh)
	    {
	      int up=loceo_neighup[eo][ieo][mu];
	      int dw=loceo_neighdw[eo][ieo][mu];
	      color_put_to_zero(out[eo][ieo]);
	      if(side==BOTH or side==UP) su3_summ_the_prod_color(out[eo][ieo],conf[eo][ieo][mu],in[!eo][up]);
	      if(side==BOTH or side==DW) su3_dag_summ_the_prod_color(out[eo][ieo],conf[!eo][dw][mu],in[!eo][dw]);
	      
	      if(side==BOTH) color_prod_double(out[eo][ieo],out[eo][ieo],0.5);
	    }
	  set_borders_invalid(out[eo]);
	}
    }
    
    //apply a single shift
    void apply_covariant_shift(color **out,quad_su3 **conf,int mu,color **in,shift_orie_t side)
    {
      for(int eo=0;eo<2;eo++) vector_reset(out[eo]);
      summ_covariant_shift(out,conf,mu,in,side);
    }
    
    //apply two-links laplacian, as per 1302.5246
    THREADABLE_FUNCTION_3ARG(two_links_laplacian, color**,out, quad_su3**,conf, color**,in)
    {
      crash("NOT USABLE YET");
      //allocate temp
      color *temp[2];
      for(int par=0;par<2;par++)
	temp[par]=nissa_malloc("temp",loc_volh+bord_volh,color);
      
      //reset out
      for(int par=0;par<2;par++) vector_reset(out[par]);
      
      //shift and summ in both orientation
      shift_orie_t orie_list[2]={UP,DW};
      for(int iorie=0;iorie<2;iorie++)
	for(int mu=0;mu<NDIM;mu++)
	  {
	    apply_covariant_shift(temp,conf,mu,in,orie_list[iorie]);
	    summ_covariant_shift(out,conf,mu,temp,orie_list[iorie]);
	  }
      
      //summ diagonal and normalize
      for(int par=0;par<2;par++)
	{
	  double_vector_summassign_double_vector_prod_double((double*)(out[par]),(double*)(in[par]),-2*NDIM,loc_volh);
	  double_vector_prodassign_double((double*)(out[par]),0.5,loc_volh);
	}
      
      //free
      for(int icopy=0;icopy<2;icopy++)
	for(int par=0;par<2;par++)
	  nissa_free(temp[icopy][par]);
    }
    THREADABLE_FUNCTION_END
    
    //apply the operator
    void apply_op_single_perm(color **out,color **temp,quad_su3 **conf,std::vector<int> &list_dir,color **in)
    {
      //make a temporary copy
      for(int eo=0;eo<2;eo++) vector_copy(temp[eo],in[eo]);
      
      for(std::vector<int>::iterator mu_it=list_dir.begin();mu_it!=list_dir.end();mu_it++)
	{
	  //write comment, copy and communicate
	  verbosity_lv2_master_printf(" shift %d\n",*mu_it);
	  if(mu_it!=list_dir.begin())
	    for(int eo=0;eo<2;eo++)
	      vector_copy(temp[eo],out[eo]);
	  
	  //make the shift
	  apply_covariant_shift(out,conf,*mu_it,temp,BOTH);
	}
    }
    
    //apply the operator summing all permutations
    void apply_op(color **out,color **single_perm,color **internal_temp,quad_su3 **conf,quad_u1 **u1b,int shift,color **in)
    {
      //make a list that can be easily permuted
      std::vector<int> list_dir;
      for(int mu=0;mu<NDIM;mu++)
	if((shift>>mu)&0x1)
	  list_dir.push_back(mu);
      std::sort(list_dir.begin(),list_dir.end());
      
      if(list_dir.size())
	{
	  add_backfield_without_stagphases_to_conf(conf,u1b);
	  
	  //summ all perms
	  int nperm=0;
	  for(int eo=0;eo<2;eo++) vector_reset(out[eo]);
	  do
	    {
	      //incrementing the number of permutations
	      verbosity_lv2_master_printf("Considering permutation %d:",nperm);
	      for(std::vector<int>::iterator it=list_dir.begin();it!=list_dir.end();it++) verbosity_lv2_master_printf(" %d",*it);
	      verbosity_lv2_master_printf("\n");
	      nperm++;
	      
	      //apply and summ
	      apply_op_single_perm(single_perm,internal_temp,conf,list_dir,in);
	      for(int eo=0;eo<2;eo++) double_vector_summassign((double*)(out[eo]),(double*)(single_perm[eo]),loc_volh*sizeof(color)/sizeof(double));
	    }
	  while(std::next_permutation(list_dir.begin(),list_dir.end()));
	  
	  //final normalization
	  for(int eo=0;eo<2;eo++) double_vector_prod_double((double*)(out[eo]),(double*)(out[eo]),1.0/nperm,loc_volh*sizeof(color)/sizeof(double));
	  
	  rem_backfield_without_stagphases_from_conf(conf,u1b);
	}
      else for(int eo=0;eo<2;eo++) vector_copy(out[eo],in[eo]);
    }
    
    //add the phases
    THREADABLE_FUNCTION_2ARG(put_stag_phases, color**,source, int,mask)
    {
      GET_THREAD_ID();
      
      //put the phases
      for(int eo=0;eo<2;eo++)
	{
	  NISSA_PARALLEL_LOOP(ieo,0,loc_volh)
	    {
	      int sign=1,ivol=loclx_of_loceo[eo][ieo];
	      for(int mu=0;mu<NDIM;mu++) sign*=1-2*(((mask>>mu)&0x1) and (glb_coord_of_loclx[ivol][mu]&0x1));
	      if(abs(sign)!=1) crash("unexpected sign %d",sign);
	      color_prod_double(source[eo][ieo],source[eo][ieo],sign);
	    }
	  set_borders_invalid(source[eo]);
	}
    }
    THREADABLE_FUNCTION_END
    
    //multiply by M^-1
    THREADABLE_FUNCTION_6ARG(mult_Minv, color**,prop, quad_su3**,conf, quad_u1**,u1b, double,m, double,residue, color**,source)
    {
      add_backfield_with_stagphases_to_conf(conf,u1b);
      inv_stD_cg(prop,conf,m,100000,residue,source);
      rem_backfield_with_stagphases_from_conf(conf,u1b);
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
    void fill_source(color **src,int twall)
    {generate_fully_undiluted_eo_source(src,RND_GAUSS,twall);}
    
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
      
      add_backfield_with_stagphases_to_conf(conf,theory_pars->backfield[iflav]);
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
      
      rem_backfield_with_stagphases_from_conf(conf,theory_pars->backfield[iflav]);
    }
    THREADABLE_FUNCTION_END
    
    //compute a density
    THREADABLE_FUNCTION_10ARG(summ_dens, complex*,dens, color**,quark, color**,temp0, color**,temp1, quad_su3**,conf, quad_u1**,backfield, int,shift, int,mask, color**,chi, color**,eta)
    {
      GET_THREAD_ID();
      
      apply_op(quark,temp0,temp1,conf,backfield,shift,chi);
      put_stag_phases(quark,mask);
      
      for(int eo=0;eo<2;eo++)
	NISSA_PARALLEL_LOOP(ieo,0,loc_volh)
	  {
	    int ivol=loclx_of_loceo[eo][ieo];
	    complex prod;
	    color_scalar_prod(prod,eta[eo][ieo],quark[eo][ieo]);
	    complex_summassign(dens[ivol],prod);
	  }
      THREAD_BARRIER();
    }
    THREADABLE_FUNCTION_END
    
    void insert_external_source_handle(complex out,spin1field **aux,int par,int ieo,int mu,void *pars)
    {if(aux) complex_copy(out,aux[par][ieo][mu]);else complex_put_to_real(out,1);}
    //insert an external current
    void insert_vector_vertex(color **out,quad_su3 **conf,theory_pars_t *theory_pars,int iflav,spin1field **curr,color **in,complex fact_fw,complex fact_bw,void(*get_curr)(complex out,spin1field **curr,int par,int ieo,int mu,void *pars),int t,void *pars)
    {
      GET_THREAD_ID();
      
      add_backfield_with_stagphases_to_conf(conf,theory_pars->backfield[iflav]);
      communicate_ev_and_od_quad_su3_borders(conf);
      communicate_ev_and_od_color_borders(in);
      if(curr) communicate_ev_and_od_spin1field_borders(curr);
      
      for(int par=0;par<2;par++)
	{
	  NISSA_PARALLEL_LOOP(ieo,0,loc_volh)
	    {
	      int ivol=loclx_of_loceo[par][ieo];
	      
	      color_put_to_zero(out[par][ieo]);
	      if(t<0  or  t>=glb_size[0]  or  glb_coord_of_loclx[ivol][0]==t)
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
		    color_subt_the_prod_complex(out[par][ieo],temp,cb);
		  }
	    }
	  set_borders_invalid(out[par]);
	}
      
      rem_backfield_with_stagphases_from_conf(conf,theory_pars->backfield[iflav]);
    }
  }
  
  std::string base_fermionic_meas_t::get_str(bool full)
  {
    std::ostringstream os;
    
    if(each!=def_each() or full) os<<" Each\t\t=\t"<<each<<"\n";
    if(after!=def_after() or full) os<<" After\t\t=\t"<<after<<"\n";
    if(path!=def_path() or full) os<<" Path\t\t=\t\""<<path.c_str()<<"\"\n";
    if(residue!=def_residue() or full) os<<" Residue\t=\t"<<residue<<"\n";
    if(ncopies!=def_ncopies() or full) os<<" NCopies\t=\t"<<ncopies<<"\n";
    if(itheory!=def_itheory() or full) os<<" ITheory\t=\t"<<itheory<<"\n";
    if(nhits!=def_nhits() or full) os<<" NHits\t\t=\t"<<nhits<<"\n";
    
    return os.str();
  }
}
