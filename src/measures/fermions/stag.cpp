#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "base/random.hpp"
#include "base/vectors.hpp"
#include "communicate/borders.hpp"
#include "geometry/geometry_eo.hpp"
#include "hmc/backfield.hpp"
#include "hmc/theory_pars.hpp"
#include "inverters/staggered/cg_invert_stD.hpp"
#include "linalgs/linalgs.hpp"
#include "linalgs/reduce.hpp"
#include "new_types/su3.hpp"
#include "measures/fermions/fermionic_meas.hpp"

#include "stag.hpp"

namespace nissa
{
  namespace stag
  {
    //summ a single shift
    void summ_covariant_shift(eo_ptr<color> out,eo_ptr<quad_su3> conf,int mu,eo_ptr<color> in,shift_orie_t side)
    {
      
      if(in==out) crash("in==out");
      
      communicate_ev_and_od_color_borders(in);
      communicate_ev_and_od_quad_su3_borders(conf);
      
      for(int eo=0;eo<2;eo++)
	{
	  NISSA_PARALLEL_LOOP(ieo,0,locVolh)
	    {
	      int up=loceo_neighup[eo][ieo][mu];
	      int dw=loceo_neighdw[eo][ieo][mu];
	      color_put_to_zero(out[eo][ieo]);
	      if(side==BOTH or side==UP) su3_summ_the_prod_color(out[eo][ieo],conf[eo][ieo][mu],in[!eo][up]);
	      if(side==BOTH or side==DW) su3_dag_summ_the_prod_color(out[eo][ieo],conf[!eo][dw][mu],in[!eo][dw]);
	      
	      if(side==BOTH) color_prod_double(out[eo][ieo],out[eo][ieo],0.5);
	    }
	  NISSA_PARALLEL_LOOP_END;
	  set_borders_invalid(out[eo]);
	}
    }
    
    //apply a single shift
    void apply_covariant_shift(eo_ptr<color> out,eo_ptr<quad_su3> conf,int mu,eo_ptr<color> in,shift_orie_t side)
    {
      for(int eo=0;eo<2;eo++) vector_reset(out[eo]);
      summ_covariant_shift(out,conf,mu,in,side);
    }
    
    //apply two-links laplacian, as per 1302.5246
    void two_links_laplacian(eo_ptr<color> out,eo_ptr<quad_su3> conf,eo_ptr<color> in)
    {
      crash("NOT USABLE YET");
      //allocate temp
      eo_ptr<color> temp;
      for(int par=0;par<2;par++)
	temp[par]=nissa_malloc("temp",locVolh+bord_volh,color);
      
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
	  double_vector_summassign_double_vector_prod_double((double*)(out[par]),(double*)(in[par]),-2*NDIM,locVolh);
	  double_vector_prodassign_double((double*)(out[par]),0.5,locVolh);
	}
      
      //free
	for(int par=0;par<2;par++)
	  nissa_free(temp[par]);
    }
    
    //apply the operator
    void apply_shift_op_single_perm(eo_ptr<color> out,eo_ptr<color> temp,eo_ptr<quad_su3> conf,std::vector<int> &list_dir,eo_ptr<color> in)
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
    void apply_shift_op(eo_ptr<color> out,eo_ptr<color> single_perm,eo_ptr<color> internal_temp,eo_ptr<quad_su3> conf,eo_ptr<quad_u1> u1b,int shift,eo_ptr<color> in)
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
	      apply_shift_op_single_perm(single_perm,internal_temp,conf,list_dir,in);
	      for(int eo=0;eo<2;eo++) double_vector_summassign((double*)(out[eo]),(double*)(single_perm[eo]),locVolh*sizeof(color)/sizeof(double));
	    }
	  while(std::next_permutation(list_dir.begin(),list_dir.end()));
	  
	  //final normalization
	  for(int eo=0;eo<2;eo++) double_vector_prod_double((double*)(out[eo]),(double*)(out[eo]),1.0/nperm,locVolh*sizeof(color)/sizeof(double));
	  
	  rem_backfield_without_stagphases_from_conf(conf,u1b);
	}
      else for(int eo=0;eo<2;eo++) vector_copy(out[eo],in[eo]);
    }
    
    //add the phases
    void put_stag_phases(eo_ptr<color> source,int mask)
    {
      
      //put the phases
      for(int eo=0;eo<2;eo++)
	{
	  NISSA_PARALLEL_LOOP(ieo,0,locVolh)
	    {
	      int sign=1,ivol=loclx_of_loceo[eo][ieo];
	      for(int mu=0;mu<NDIM;mu++) sign*=1-2*(((mask>>mu)&0x1) and (glbCoordOfLoclx[ivol][mu]&0x1));
	      color_prod_double(source[eo][ieo],source[eo][ieo],sign);
	    }
	  NISSA_PARALLEL_LOOP_END;
	  set_borders_invalid(source[eo]);
	}
    }
    
    //multiply by M^-1
    void mult_Minv(eo_ptr<color> prop,eo_ptr<quad_su3> conf,eo_ptr<quad_u1> u1b,double m,double residue,eo_ptr<color> source)
    {
      add_backfield_with_stagphases_to_conf(conf,u1b);
      inv_stD_cg(prop,conf,m,100000,residue,source);
      rem_backfield_with_stagphases_from_conf(conf,u1b);
    }
    void mult_Minv(eo_ptr<color> prop,eo_ptr<quad_su3> conf,theory_pars_t *pars,int iflav,double residue,eo_ptr<color> source)
    {mult_Minv(prop,conf,pars->backfield[iflav],pars->quarks[iflav].mass,residue,source);}
    
    //compute the matrix element of the derivative of the dirac operator between two vectors
    //forward and backward derivative are stored separately, for a reason
    void compute_fw_bw_der_mel(complex *res_fw_bw,eo_ptr<color> left,eo_ptr<quad_su3> conf,int mu,eo_ptr<color> right,complex *point_result)
    {
      
      communicate_ev_and_od_color_borders(left);
      communicate_ev_and_od_quad_su3_borders(conf);
      communicate_ev_and_od_color_borders(right);
      
      for(int fw_bw=0;fw_bw<2;fw_bw++)
	{
	  vector_reset(point_result);
	  for(int par=0;par<2;par++)
	    NISSA_PARALLEL_LOOP(ieo,0,locVolh)
	      {
		eo_ptr<color> right_fw_bw[2]={right,left};
		
		color v;
		unsafe_su3_prod_color(v,conf[par][ieo][mu],right_fw_bw[fw_bw][!par][loceo_neighup[par][ieo][mu]]);
		complex t;
		if(fw_bw==0) color_scalar_prod(t,right_fw_bw[!fw_bw][par][ieo],v);
		else         color_scalar_prod(t,v,right_fw_bw[!fw_bw][par][ieo]);
		complex_summassign(point_result[loclx_of_loceo[par][ieo]],t);
	      }
	  NISSA_PARALLEL_LOOP_END;
	  THREAD_BARRIER();
	  glb_reduce(&res_fw_bw[fw_bw],point_result,locVol);
	  
	  //DEB_STAG("fw_bw=%d mu=%d, RE=%lg IM=%lg\n",fw_bw,mu,res_fw_bw[fw_bw][RE],res_fw_bw[fw_bw][IM]);
	  
	}
    }
    
    //fill a source
    void fill_source(eo_ptr<color> src,int twall,rnd_t noise_type)
    {
      generate_fully_undiluted_eo_source(src,noise_type,twall);
    }
    
    //take the trace between A^dag and B
    void summ_the_trace(double* out,complex* point_result,eo_ptr<color>  A,eo_ptr<color>  B)
    {
      
      //compute results for single points
      vector_reset(point_result);
      for(int par=0;par<2;par++)
	NISSA_PARALLEL_LOOP(ieo,0,locVolh)
	  for(int ic=0;ic<3;ic++)
	    complex_summ_the_conj1_prod(point_result[loclx_of_loceo[par][ieo]],A[par][ieo][ic],B[par][ieo][ic]);
      NISSA_PARALLEL_LOOP_END;
      THREAD_BARRIER();
      
      //final reduction
      complex temp;
      glb_reduce(&temp,point_result,locVol);
      if(IS_MASTER_THREAD) complex_summassign(out,temp);
    }
    
    //multiply by the derivative of M w.r.t mu
    void mult_dMdmu(eo_ptr<color> out,theory_pars_t* theory_pars,eo_ptr<quad_su3> conf,int iflav,int ord,eo_ptr<color> in)
    {
      
      if(ord==0) crash("makes no sense to call with order zero");
      
      add_backfield_with_stagphases_to_conf(conf,theory_pars->backfield[iflav]);
      communicate_ev_and_od_quad_su3_borders(conf);
      communicate_ev_and_od_color_borders(in);
      
      for(int par=0;par<2;par++)
	{
	  NISSA_PARALLEL_LOOP(ieo,0,locVolh)
	    {
	      color temp;
	      unsafe_su3_prod_color(temp,conf[par][ieo][0],in[!par][loceo_neighup[par][ieo][0]]);
	      int idw=loceo_neighdw[par][ieo][0];
	      if(ord%2==0) su3_dag_subt_the_prod_color(temp,conf[!par][idw][0],in[!par][idw]);
	      else         su3_dag_summ_the_prod_color(temp,conf[!par][idw][0],in[!par][idw]);
	      color_prod_double(out[par][ieo],temp,0.5);
	    }
	  NISSA_PARALLEL_LOOP_END;
	  set_borders_invalid(out[par]);
	}
      
      rem_backfield_with_stagphases_from_conf(conf,theory_pars->backfield[iflav]);
    }
    
    //compute a density
    void summ_dens(complex* dens,eo_ptr<color> quark,eo_ptr<color> temp0,eo_ptr<color> temp1,eo_ptr<quad_su3> conf,eo_ptr<quad_u1> backfield,int shift,int mask,eo_ptr<color> chi,eo_ptr<color> eta)
    {
      
      apply_shift_op(quark,temp0,temp1,conf,backfield,shift,chi);
      put_stag_phases(quark,mask);
      
      for(int eo=0;eo<2;eo++)
	NISSA_PARALLEL_LOOP(ieo,0,locVolh)
	  {
	    int ivol=loclx_of_loceo[eo][ieo];
	    complex prod;
	    color_scalar_prod(prod,eta[eo][ieo],quark[eo][ieo]);
	    complex_summassign(dens[ivol],prod);
	  }
      NISSA_PARALLEL_LOOP_END;
      THREAD_BARRIER();
    }
    
    void insert_external_source_handle(complex out,eo_ptr<spin1field> aux,int par,int ieo,int mu,void *pars)
    {if(aux[0]) complex_copy(out,aux[par][ieo][mu]);else complex_put_to_real(out,1);}
    //insert an external current
    void insert_vector_vertex(eo_ptr<color> out,eo_ptr<quad_su3> conf,theory_pars_t *theory_pars,int iflav,eo_ptr<spin1field> curr,eo_ptr<color> in,complex fact_fw,complex fact_bw,void(*get_curr)(complex out,eo_ptr<spin1field> curr,int par,int ieo,int mu,void *pars),int t,void *pars)
    {
      
      add_backfield_with_stagphases_to_conf(conf,theory_pars->backfield[iflav]);
      communicate_ev_and_od_quad_su3_borders(conf);
      communicate_ev_and_od_color_borders(in);
      if(curr[0]) communicate_ev_and_od_spin1field_borders(curr);
      
      for(int par=0;par<2;par++)
	{
	  NISSA_PARALLEL_LOOP(ieo,0,locVolh)
	    {
	      int ivol=loclx_of_loceo[par][ieo];
	      
	      color_put_to_zero(out[par][ieo]);
	      if(t<0  or  t>=glbSize[0]  or  glbCoordOfLoclx[ivol][0]==t)
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
	  NISSA_PARALLEL_LOOP_END;
	  set_borders_invalid(out[par]);
	}
      
      rem_backfield_with_stagphases_from_conf(conf,theory_pars->backfield[iflav]);
    }
  }
}
