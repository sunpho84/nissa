#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "base/vectors.hpp"
#include "geometry/geometry_eo.hpp"
#include "new_types/su3_op.hpp"
#include "linalgs/linalgs.hpp"
#include "threads/threads.hpp"
#include "operations/su3_paths/clover_term.hpp"
#include "hmc/quark_pars.hpp"

namespace nissa
{
  void clover_det_force(eo_ptr<quad_su3> F,std::vector<quark_content_t> quark_content,eo_ptr<quad_su3> eo_conf)
  {
    bool need=false;
    for(auto& q : quark_content)
      need|=(q.cSW!=0);
    
    if(need)
      {
	
	//Prepare clover
	eo_ptr<clover_term_t> Cl;
	for(int eo=0;eo<2;eo++)
	  Cl[eo]=nissa_malloc("Cl",locVolh.nastyConvert(),clover_term_t);
	chromo_operator(Cl,eo_conf);
	
	eo_ptr<inv_clover_term_t> invCl;
	for(int eo=0;eo<2;eo++)
	  invCl[eo]=nissa_malloc("invCl",locVolh.nastyConvert(),inv_clover_term_t);
	
	as2t_su3 *insertion=nissa_malloc("insertion",locVolhWithBordAndEdge.nastyConvert(),as2t_su3);
	for(auto& q : quark_content)
	  {
	    const double cSW=q.cSW;
	    
	    if(cSW)
	      {
		chromo_operator_include_cSW(Cl,q.cSW);
		
		for(int eo=0;eo<2;eo++)
		  invert_twisted_clover_term(invCl[eo],q.mass,q.kappa,Cl[eo]);
		
		/////////////////////////////////////////////////////////////////
		
		NISSA_PARALLEL_LOOP(jeo,0,locVolh)
		  {
		    FOR_ALL_DIRECTIONS(mu)
		      for(Direction nu=mu+1;nu<NDIM;nu++)
			{
			  int ipair=edge_numb[mu.nastyConvert()][nu.nastyConvert()];
			  dirac_matr m=dirac_prod(base_gamma[igamma_of_mu(mu)],base_gamma[igamma_of_mu(nu)]);
			  
			  su3& ins=insertion[jeo.nastyConvert()][ipair];
			  
			  for(int ic1=0;ic1<NCOL;ic1++)
			    for(int ic2=0;ic2<NCOL;ic2++)
			      {
				complex_put_to_zero(ins[ic1][ic2]);
				
				for(int x_high_low=0;x_high_low<2;x_high_low++)
				  for(int iw=0;iw<NDIRAC/2;iw++)
				    {
				      int id=2*x_high_low+iw;
				      complex& c=m.entr[id];
				      int jd=m.pos[id];
				      int jw=jd-2*x_high_low;
				      
				      complex_summ_the_prod(ins[ic1][ic2],c,invCl[EVN][jeo.nastyConvert()][x_high_low][jw][ic1][iw][ic2]);
				    }
			      }
			  
			  su3_anti_hermitian_part(ins,ins);
			}
		  }
		NISSA_PARALLEL_LOOP_END;
		
		set_borders_invalid(insertion);
		
		//communicate_e o_as2t_su3_edges(cl_insertion[]);
		
		FOR_BOTH_PARITIES(par)
		  NISSA_PARALLEL_LOOP(ieo,0,locVolh)
		    {
		      FOR_ALL_DIRECTIONS(mu)
			{
			  su3 contr;
			  su3_put_to_zero(contr);
			  
			  for(int inu=0;inu<NDIM-1;inu++)
			    {
			      const Direction nu=perp_dir[mu.nastyConvert()][inu];
			      
			      const LocEoSite xpmu=loceo_neighup(par,ieo,mu);
			      const LocEoSite xmnu=loceo_neighdw(par,ieo,nu);
			      const LocEoSite xpnu=loceo_neighup(par,ieo,nu);
			      const LocEoSite xpmumnu=loceo_neighdw((1-par),xpmu,nu);
			      const LocEoSite xpmupnu=loceo_neighup((1-par),xpmu,nu);
			      
			      int ipair=edge_numb[mu.nastyConvert()][nu.nastyConvert()];
			      
			      for(int i=0;i<2;i++)
				{
				  su3 u;
				  
				  double sign;
				  if(mu<nu) sign=+1.0;
				  else       sign=-1.0;
				  
				  su3_put_to_diag(u,sign);
				  if(i==0 and par==ODD) safe_su3_prod_su3(u,u,insertion[xpmu.nastyConvert()][ipair]);
				  safe_su3_prod_su3(u,u,eo_conf[(1-par).nastyConvert()][xpmu.nastyConvert()][nu.nastyConvert()]);
				  if(i==0 and par==EVN) safe_su3_prod_su3(u,u,insertion[xpmupnu.nastyConvert()][ipair]);
				  safe_su3_prod_su3_dag(u,u,eo_conf[(1-par).nastyConvert()][xpnu.nastyConvert()][mu.nastyConvert()]);
				  if(i==1 and par==ODD) safe_su3_prod_su3(u,u,insertion[xpnu.nastyConvert()][ipair]);
				  safe_su3_prod_su3_dag(u,u,eo_conf[par][ieo.nastyConvert()][nu.nastyConvert()]);
				  if(i==1 and par==EVN) safe_su3_prod_su3(u,u,insertion[ieo.nastyConvert()][ipair]);
				  
				  su3_summassign(contr,u);
				  
				  su3 v;
				  
				  su3_put_to_diag(v,sign);
				  if(i==0 and par==ODD) safe_su3_prod_su3(v,v,insertion[xpmu.nastyConvert()][ipair]);
				  safe_su3_prod_su3_dag(v,v,eo_conf[par][xpmumnu.nastyConvert()][nu.nastyConvert()]);
				  if(i==0 and par==EVN) safe_su3_prod_su3(v,v,insertion[xpmumnu.nastyConvert()][ipair]);
				  safe_su3_prod_su3_dag(v,v,eo_conf[(1-par).nastyConvert()][xmnu.nastyConvert()][mu.nastyConvert()]);
				  if(i==1 and par==ODD) safe_su3_prod_su3(v,v,insertion[xmnu.nastyConvert()][ipair]);
				  safe_su3_prod_su3(v,v,eo_conf[(1-par).nastyConvert()][xmnu.nastyConvert()][nu.nastyConvert()]);
				  if(i==1 and par==EVN) safe_su3_prod_su3(v,v,insertion[ieo.nastyConvert()][ipair]);
				  
				  su3_subtassign(contr,v);
				}
			    }
			  
			  su3_prodassign_double(contr,-cSW/4);
			}
		    }
		NISSA_PARALLEL_LOOP_END;
		
		chromo_operator_remove_cSW(Cl,q.cSW);
	      }
	  }
	nissa_free(insertion);
	
	for(int eo=0;eo<2;eo++)
	  {
	    nissa_free(Cl[eo]);
	    nissa_free(invCl[eo]);
	  }
      }
  }
}
