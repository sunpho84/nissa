#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "base/bench.hpp"
#include "base/vectors.hpp"
#include "communicate/borders.hpp"
#include "communicate/edges.hpp"
#include "dirac_operators/tmD_eoprec/dirac_operator_tmD_eoprec.hpp"
#include "dirac_operators/tmclovD_eoprec/dirac_operator_tmclovD_eoprec.hpp"
#include "hmc/backfield.hpp"
#include "inverters/twisted_clover/cgm_invert_tmclovDkern_eoprec_square_portable.hpp"
#include "new_types/su3_op.hpp"
#include "threads/threads.hpp"

namespace nissa
{
  /// Derivative of xQy
  void get_point_twisted_force(su3 out,spincolor *a[2],spincolor *b[2],int eo,int ieo,int dir)
  {
    int ineoup=loceo_neighup[eo][ieo][dir];
    
    spincolor temp;
    spincolor_copy(temp,a[!eo][ineoup]);
    dirac_subt_the_prod_spincolor(temp,base_gamma+igamma_of_mu[dir],a[!eo][ineoup]);
    safe_dirac_prod_spincolor(temp,base_gamma+5,temp);
    
    for(int ic1=0;ic1<NCOL;ic1++)
      for(int ic2=0;ic2<NCOL;ic2++)
	{
	  complex& o=out[ic1][ic2];
	  complex_put_to_zero(o);
	  
	  for(int id=0;id<NDIRAC;id++)
	    complex_subt_the_conj2_prod(o,temp[id][ic1],b[eo][ieo][id][ic2]);
	}
  }
  
  // implement appendix B of https://arxiv.org/pdf/0905.3331.pdf
  THREADABLE_FUNCTION_10ARG(summ_the_roottm_clov_eoimpr_quark_force, quad_su3**,F, quad_su3**,eo_conf, double,kappa, double,cSW, clover_term_t*,Cl_odd, inv_clover_term_t*,invCl_evn, spincolor*,phi_o, quad_u1**,u1b, rat_approx_t*,appr, double,residue)
  {
    GET_THREAD_ID();
    
    START_TIMING(quark_force_over_time,nquark_force_over);
    
    //allocate each terms of the expansion
    spincolor *Y[2][appr->degree()],*X[2][appr->degree()],*temp=nissa_malloc("temp",loc_volh+bord_volh,spincolor);
    for(int eo=0;eo<2;eo++)
      for(int iterm=0;iterm<appr->degree();iterm++)
	{
	  Y[eo][iterm]=nissa_malloc("Y",loc_volh+bord_volh,spincolor);
	  X[eo][iterm]=nissa_malloc("X",loc_volh+bord_volh,spincolor);
	}
    
    //add the background fields
    add_backfield_without_stagphases_to_conf(eo_conf,u1b);
    
    //invert the various terms
    STOP_TIMING(quark_force_over_time);
    // eq. B.8a
    inv_tmclovDkern_eoprec_square_portable_run_hm_up_to_comm_prec(X[ODD],eo_conf,kappa,Cl_odd,invCl_evn,appr->poles.data(),appr->degree(),10000000,residue,phi_o);
    UNPAUSE_TIMING(quark_force_over_time);
    
    ////////////////////
    
    for(int iterm=0;iterm<appr->degree();iterm++)
      {
	//eq. B.8b: Y_o = Qhat^- X_o
	tmclovDkern_eoprec_eos(Y[ODD][iterm],temp,eo_conf,kappa,Cl_odd,invCl_evn,true,appr->poles[iterm],X[ODD][iterm]);
	
	tmn2Deo_eos(temp,eo_conf,X[ODD][iterm]); // temp = - 2 * D_eo * X_o
	inv_tmclovDee_or_oo_eos(X[EVN][iterm],invCl_evn,true,temp); // X_e = M_ee+^-1 * temp = - 2 * M_ee+^-1 * D_eo * X_o
	double_vector_prodassign_double((double*)(X[EVN][iterm]),0.5,loc_volh*sizeof(spincolor)/sizeof(double)); // X_e = 0.5 * X_e = - M_ee+^-1 * D_eo * X_o
	
	tmn2Deo_eos(temp,eo_conf,Y[ODD][iterm]); // temp = - 2 * D_eo * Y_o
	inv_tmclovDee_or_oo_eos(Y[EVN][iterm],invCl_evn,false,temp); // Y_e = M_ee-^-1 * temp = - 2 * M_ee-^-1 * D_eo * Y_o
	double_vector_prodassign_double((double*)(Y[EVN][iterm]),0.5,loc_volh*sizeof(spincolor)/sizeof(double)); // Y_e = 0.5 * Y_e = - M_ee-^-1 * D_eo * Y_o
      }
    
    //communicate borders (could be improved...)
    for(int eo=0;eo<2;eo++)
      for(int iterm=0;iterm<appr->degree();iterm++)
	for(auto v : {X[eo][iterm],Y[eo][iterm]})
	  communicate_ev_or_od_spincolor_borders(v,eo);
    
    //conclude the calculation of the fermionic force
    for(int iterm=0;iterm<appr->degree();iterm++)
      for(int eo=0;eo<2;eo++)
	NISSA_PARALLEL_LOOP(ieo,0,loc_volh)
	  for(int mu=0;mu<NDIM;mu++)
	    {
	      su3& out=F[eo][ieo][mu];
	      
	      int ineoup=loceo_neighup[eo][ieo][mu];
	      
	      for(int i=0;i<2;i++)
		{
		  spincolor& a=((i==0)?X[!eo]:Y[!eo])[iterm][ineoup];
		  spincolor& b=((i==0)?Y[eo]:X[eo])[iterm][ieo];
		  
		  spincolor temp;
		  spincolor_copy(temp,a);
		  dirac_subt_the_prod_spincolor(temp,base_gamma+igamma_of_mu[mu],a);
		  safe_dirac_prod_spincolor(temp,base_gamma+5,temp);
		  
		  for(int ic1=0;ic1<NCOL;ic1++)
		    for(int ic2=0;ic2<NCOL;ic2++)
		      for(int id=0;id<NDIRAC;id++)
			complex_subt_the_conj2_prod(out[ic1][ic2],temp[id][ic1],b[id][ic2]);
		}
	    }
    NISSA_PARALLEL_LOOP_END;
    
    if(cSW!=0)
      {
	as2t_su3 *cl_insertion[2];
	for(int eo=0;eo<2;eo++)
	  cl_insertion[eo]=nissa_malloc("insertion",loc_volh+bord_volh+edge_volh,as2t_su3);
	
	for(int iterm=0;iterm<appr->degree();iterm++)
	  for(int eo=0;eo<2;eo++)
	    {
	      NISSA_PARALLEL_LOOP(jeo,0,loc_volh)
		{
		  for(int mu=0;mu<NDIM;mu++)
		    for(int nu=mu+1;nu<NDIM;nu++)
		      {
			int ipair=edge_numb[mu][nu];
			dirac_matr m=dirac_prod(base_gamma[5],dirac_prod(base_gamma[igamma_of_mu[mu]],base_gamma[igamma_of_mu[nu]]));
			
			su3& ins=cl_insertion[eo][jeo][ipair];
			spincolor tempX,tempY;
			unsafe_dirac_prod_spincolor(tempX,&m,X[eo][iterm][jeo]);
			unsafe_dirac_prod_spincolor(tempY,&m,Y[eo][iterm][jeo]);
			
			su3_put_to_zero(ins);
			
			for(int ic1=0;ic1<NCOL;ic1++)
			  for(int ic2=0;ic2<NCOL;ic2++)
			    for(int id=0;id<NDIRAC;id++)
			      {
				complex_summ_the_conj2_prod(ins[ic1][ic2],tempY[id][ic1],X[eo][iterm][jeo][id][ic2]);
				complex_summ_the_conj2_prod(ins[ic1][ic2],tempX[id][ic1],Y[eo][iterm][jeo][id][ic2]);
			      }
		      }
		}
	      NISSA_PARALLEL_LOOP_END;
	      set_borders_invalid(cl_insertion[eo]);
	    }
	
	communicate_eo_as2t_su3_edges(cl_insertion);
	
	for(int eo=0;eo<2;eo++)
	  NISSA_PARALLEL_LOOP(ieo,0,loc_volh)
	    for(int dir=0;dir<NDIM;dir++)
	      {
		su3 stap;
		su3_put_to_zero(stap);
		
		for(int inu=0;inu<NDIM-1;inu++)
		  {
		    int nu=perp_dir[dir][inu];
		    
		    int xpmu=loceo_neighup[eo][ieo][dir];
		    int xmnu=loceo_neighdw[eo][ieo][nu];
		    int xpnu=loceo_neighup[eo][ieo][nu];
		    int xpmumnu=loceo_neighdw[!eo][xpmu][nu];
		    int xpmupnu=loceo_neighup[!eo][xpmu][nu];
		    
		    int ipair=edge_numb[dir][nu];
		    
		    for(int i=0;i<4;i++)
		      {
			double sign;
			if(dir<nu) sign=+1.0;
			else       sign=-1.0;
			
			su3 u;
			
			su3_put_to_diag(u,sign);
			if(i==0) safe_su3_prod_su3(u,u,cl_insertion[!eo][xpmu][ipair]);
			safe_su3_prod_su3(u,u,eo_conf[!eo][xpmu][nu]);
			if(i==1) safe_su3_prod_su3(u,u,cl_insertion[eo][xpmupnu][ipair]);
			safe_su3_prod_su3_dag(u,u,eo_conf[!eo][xpnu][dir]);
			if(i==2) safe_su3_prod_su3(u,u,cl_insertion[!eo][xpnu][ipair]);
			safe_su3_prod_su3_dag(u,u,eo_conf[eo][ieo][nu]);
			if(i==3) safe_su3_prod_su3(u,u,cl_insertion[eo][ieo][ipair]);
			su3_summassign(stap,u);
			
			su3 v;
			
			su3_put_to_diag(v,sign);
			if(i==0) safe_su3_prod_su3(v,v,cl_insertion[!eo][xpmu][ipair]);
			safe_su3_prod_su3_dag(v,v,eo_conf[eo][xpmumnu][nu]);
			if(i==1) safe_su3_prod_su3(v,v,cl_insertion[eo][xpmumnu][ipair]);
			safe_su3_prod_su3_dag(v,v,eo_conf[!eo][xmnu][dir]);
			if(i==2) safe_su3_prod_su3(v,v,cl_insertion[!eo][xmnu][ipair]);
			safe_su3_prod_su3(v,v,eo_conf[!eo][xmnu][nu]);
			if(i==3) safe_su3_prod_su3(v,v,cl_insertion[eo][ieo][ipair]);
			su3_subtassign(stap,v);
		      }
		    
		    su3_summ_the_prod_double(F[eo][ieo][dir],stap,-cSW/8);
		  }
	      }
	NISSA_PARALLEL_LOOP_END;
	
	for(int eo=0;eo<2;eo++)
	  nissa_free(cl_insertion[eo]);
      }
    #warning extra fase u1?
    //remove the background fields
    rem_backfield_without_stagphases_from_conf(eo_conf,u1b);
    
    //free
    for(int eo=0;eo<2;eo++)
      for(int iterm=0;iterm<appr->degree();iterm++)
	{
	  nissa_free(X[eo][iterm]);
	  nissa_free(Y[eo][iterm]);
	}
    nissa_free(temp);
    
    STOP_TIMING(quark_force_over_time);
  }
  THREADABLE_FUNCTION_END
}
