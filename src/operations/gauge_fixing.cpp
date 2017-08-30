#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <math.h>
#include <string.h>

#include "base/debug.hpp"
#include "base/random.hpp"
#include "base/thread_macros.hpp"
#include "base/vectors.hpp"
#include "communicate/borders.hpp"
#include "geometry/geometry_eo.hpp"
#include "geometry/geometry_lx.hpp"
#include "linalgs/linalgs.hpp"
#include "new_types/complex.hpp"
#include "new_types/su3_op.hpp"
#include "operations/fft.hpp"
#include "routines/ios.hpp"
#include "routines/mpi_routines.hpp"

#ifdef USE_THREADS
 #include "routines/thread.hpp"
#endif

namespace nissa
{
  //apply the passed transformation to the point
  void local_gauge_transform(quad_su3 *conf,su3 g,int ivol)
  {
    // for each dir...
    for(int mu=0;mu<NDIM;mu++)
      {
        int b=loclx_neighdw[ivol][mu];
        
        //perform local gauge transform
        safe_su3_prod_su3(conf[ivol][mu],g,conf[ivol][mu]);
        safe_su3_prod_su3_dag(conf[b][mu],conf[b][mu],g);
      }
  }
  
  //apply a gauge transformation to the conf
  THREADABLE_FUNCTION_3ARG(gauge_transform_conf, quad_su3*,uout, su3*,g, quad_su3*,uin)
  {
    GET_THREAD_ID();
    
    //communicate borders
    communicate_lx_su3_borders(g);
    
    //transform
    su3 temp;
    NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
      for(int mu=0;mu<NDIM;mu++)
	{
	  unsafe_su3_prod_su3_dag(temp,uin[ivol][mu],g[loclx_neighup[ivol][mu]]);
	  unsafe_su3_prod_su3(uout[ivol][mu],g[ivol],temp);
	}
    
    //invalidate borders
    set_borders_invalid(uout);
  }
  THREADABLE_FUNCTION_END
  //e/o version
  THREADABLE_FUNCTION_3ARG(gauge_transform_conf, quad_su3**,uout, su3**,g, quad_su3**,uin)
  {
    GET_THREAD_ID();
    
    //communicate borders
    communicate_ev_and_od_su3_borders(g);
    
    //transform
    su3 temp;
    for(int par=0;par<2;par++)
      NISSA_PARALLEL_LOOP(ivol,0,loc_volh)
	for(int mu=0;mu<NDIM;mu++)
	  {
	    unsafe_su3_prod_su3_dag(temp,uin[par][ivol][mu],g[!par][loceo_neighup[par][ivol][mu]]);
	    unsafe_su3_prod_su3(uout[par][ivol][mu],g[par][ivol],temp);
	  }
    
    //invalidate borders
    set_borders_invalid(uout[0]);
    set_borders_invalid(uout[1]);
  }
  THREADABLE_FUNCTION_END
  
  //transform a color field
  THREADABLE_FUNCTION_3ARG(gauge_transform_color, color**,out, su3**,g, color**,in)
  {
    GET_THREAD_ID();
    
    //communicate borders
    communicate_ev_and_od_su3_borders(g);
    
    //transform
    for(int par=0;par<2;par++)
      NISSA_PARALLEL_LOOP(ivol,0,loc_volh)
	safe_su3_prod_color(out[par][ivol],g[par][ivol],in[par][ivol]);
    
    //invalidate borders
    set_borders_invalid(out[0]);
    set_borders_invalid(out[1]);
  }
  THREADABLE_FUNCTION_END
  
  //determine the gauge transformation bringing to temporal gauge with T-1 timeslice diferent from id
  void find_temporal_gauge_fixing_matr(su3 *fixm,quad_su3 *u)
  {
    int loc_slice_area=loc_size[1]*loc_size[2]*loc_size[3];
    su3 *buf=NULL;
    
    //if the number of ranks in the 0 dir is greater than 1 allocate room for border
    if(nrank_dir[0]>1) buf=nissa_malloc("buf",loc_slice_area,su3);
    
    //if we are on first rank slice put to identity the t=0 slice, otherwise receive it from previous rank slice
    if(rank_coord[0]==0)
      {
	NISSA_LOC_VOL_LOOP(ivol)
	  if(glb_coord_of_loclx[ivol][0]==0)
	    su3_put_to_id(fixm[ivol]);
      }
    else
      if(nrank_dir[0]>1)
	MPI_Recv((void*)fixm,loc_slice_area,MPI_SU3,rank_neighdw[0],252,cart_comm,MPI_STATUS_IGNORE);
    
    //now go ahead along t
    int c[NDIM];
    //loop over spatial slice
    for(c[1]=0;c[1]<loc_size[1];c[1]++)
      for(c[2]=0;c[2]<loc_size[2];c[2]++)
	for(c[3]=0;c[3]<loc_size[3];c[3]++)
	  {
	    //bulk
	    for(c[0]=1;c[0]<loc_size[0];c[0]++)
	      {
		int icurr=loclx_of_coord(c);
		c[0]--;int iback=loclx_of_coord(c);c[0]++;
		
		unsafe_su3_prod_su3(fixm[icurr],fixm[iback],u[iback][0]);
	      }
	    //border
	    if(nrank_dir[0]>1)
	      {
		c[0]=loc_size[0]-1;int iback=loclx_of_coord(c);
		c[0]=0;int icurr=loclx_of_coord(c);
		
		unsafe_su3_prod_su3(buf[icurr],fixm[iback],u[iback][0]);
	      }
	    
	  }
    
    //if we are not on last slice of rank send g to next slice
    if(rank_coord[0]!=(nrank_dir[0]-1) && nrank_dir[0]>1)
      MPI_Send((void*)buf,loc_slice_area,MPI_SU3,rank_neighup[0],252,cart_comm);
    
    if(nrank_dir[0]>1) nissa_free(buf);
  }
  
  ////////////////////////////////////// Landau or Coulomb gauges ///////////////////////////////////////////////////////
  
  //compute the functional on a single point
  double compute_Landau_or_Coulomb_functional(quad_su3 *conf,int ivol,int start_mu)
  {
    double F=0;
    
    for(int mu=start_mu;mu<NDIM;mu++)
      {
	F-=su3_real_trace(conf[ivol][mu]);
	F-=su3_real_trace(conf[loclx_neighdw[ivol][mu]][mu]);
      }
    
    return F;
  }
  
  //derivative of the functional
  void compute_Landau_or_Coulomb_functional_der(su3 out,quad_su3 *conf,int ivol,int start_mu)
  {
    su3_put_to_zero(out);
    
    for(int mu=start_mu;mu<NDIM;mu++)
      {
	su3_summassign(out,conf[ivol][mu]);
	su3_summassign_su3_dag(out,conf[loclx_neighdw[ivol][mu]][mu]);
      }
  }
  
  //compute the functional that gets minimised
  double compute_Landau_or_Coulomb_functional(quad_su3 *conf,int start_mu)
  {
    GET_THREAD_ID();
    
    double *loc_F=nissa_malloc("loc_F",loc_vol,double);
    
    NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
      {
	loc_F[ivol]=0;
	for(int mu=start_mu;mu<NDIM;mu++)
	  loc_F[ivol]-=su3_real_trace(conf[ivol][mu]);
      }
    THREAD_BARRIER();
    
    //collapse
    double F;
    double_vector_glb_collapse(&F,loc_F,loc_vol);
    nissa_free(loc_F);
    
    return F;
  }
  
  //compute the quality of the Landau or Coulomb gauge fixing
  double compute_Landau_or_Coulomb_gauge_fixing_quality(quad_su3 *conf,int start_mu)
  {
    GET_THREAD_ID();
    
    communicate_lx_quad_su3_borders(conf);
    
    double *loc_omega=nissa_malloc("loc_omega",loc_vol,double);
    
    NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
      {
	su3 delta;
	su3_put_to_zero(delta);
	
	for(int mu=start_mu;mu<NDIM;mu++)
	  {
	    su3_subtassign(delta,conf[ivol][mu]);
	    su3_summassign(delta,conf[loclx_neighdw[ivol][mu]][mu]);
	  }
	
	//take 2 the traceless anti-hermitian part
	su3 delta_TA;
	unsafe_su3_traceless_anti_hermitian_part(delta_TA,delta);
	loc_omega[ivol]=4*su3_norm2(delta_TA);
	if(ivol==0) master_printf("c: %.16lg\n",loc_omega[ivol]);
      }
    THREAD_BARRIER();
    
    //global reduction
    double omega;
    double_vector_glb_collapse(&omega,loc_omega,loc_vol);
    nissa_free(loc_omega);
    
    return omega/glb_vol/NCOL;
  }
  
  //do all the fixing
  void Landau_or_Coulomb_gauge_fix(quad_su3 *conf,su3 *fixer,int start_mu,double over_relax_prob)
  {
    GET_THREAD_ID();
    
    //#define DEBUG
    
#ifdef DEBUG

    double eps=1e-4;
    
    quad_su3 *ori_conf=nissa_malloc("ori_conf",loc_vol+bord_vol+edge_vol,quad_su3);
    vector_copy(ori_conf,conf);
    double act_ori=compute_Landau_or_Coulomb_functional(conf,start_mu);
    
    //store derivative
    su3 nu_plus,nu_minus;
    su3_put_to_zero(nu_plus);
    su3_put_to_zero(nu_minus);
    
    for(int igen=0;igen<NCOL*NCOL-1;igen++)
      {
	//prepare increment and change
	su3 ba;
	su3_prod_double(ba,gell_mann_matr[igen],eps/2);
	su3 exp_mod;
	safe_hermitian_exact_i_exponentiate(exp_mod,ba);
	
	//change -, compute action
	vector_copy(conf,ori_conf);
	local_gauge_transform(conf,exp_mod,0);
	double act_minus=compute_Landau_or_Coulomb_functional(conf,start_mu);
	
	//change +, compute action
	vector_copy(conf,ori_conf);
	safe_su3_hermitian(exp_mod,exp_mod);
	local_gauge_transform(conf,exp_mod,0);
	double act_plus=compute_Landau_or_Coulomb_functional(conf,start_mu);
	
	//set back everything
	vector_copy(conf,ori_conf);
	
	//printf("plus: %+016.016le, ori: %+16.16le, minus: %+16.16le, eps: %lg\n",act_plus,act_ori,act_minus,eps);
	double gr_plus=-(act_plus-act_ori)/eps;
	double gr_minus=-(act_ori-act_minus)/eps;
	su3_summ_the_prod_idouble(nu_plus,gell_mann_matr[igen],gr_plus);
	su3_summ_the_prod_idouble(nu_minus,gell_mann_matr[igen],gr_minus);
      }
    
    //take the average
    su3 nu;
    su3_summ(nu,nu_plus,nu_minus);
    su3_prodassign_double(nu,0.5);
    
    su3 an_temp;
    compute_Landau_or_Coulomb_functional_der(an_temp,conf,0,start_mu);
    su3 an;
    unsafe_su3_traceless_anti_hermitian_part(an,an_temp);
    
    master_printf("Comparing\n");
    master_printf("an\n");
    su3_print(an);
    master_printf("nu+\n");
    su3_print(nu_plus);
    master_printf("nu-\n");
    su3_print(nu_minus);
    master_printf("nu\n");
    su3_print(nu);
    su3 diff;
    su3_subt(diff,an,nu_plus);
    master_printf("Norm of the difference+: %lg\n",sqrt(su3_norm2(diff)));
    su3_subt(diff,an,nu_minus);
    master_printf("Norm of the difference-: %lg\n",sqrt(su3_norm2(diff)));
    su3_subt(diff,an,nu);
    master_printf("Norm of the difference: %lg\n",sqrt(su3_norm2(diff)));
    crash("pui");
#endif
    
    
    for(int eo=0;eo<2;eo++)
      {
	NISSA_PARALLEL_LOOP(ieo,0,loc_volh)
	  {
	    int ivol=loclx_of_loceo[eo][ieo];
	    
	    //compute the derivative
	    su3 temp;
	    compute_Landau_or_Coulomb_functional_der(temp,conf,ivol,start_mu);
	    
	    //dagger
	    su3 ref;
	    unsafe_su3_hermitian(ref,temp);
	    
	    //find the link that maximizes the trace
	    su3 g;
	    su3_unitarize_maximal_trace_projecting(g,ref);
	    
	    //square probabilistically
	    double p=rnd_get_unif(loc_rnd_gen+ivol,0,1);
	    if(p<over_relax_prob) safe_su3_prod_su3(g,g,g);
	    
	    //transform
	    local_gauge_transform(conf,g,ivol);
	    
	    //store the change
	    safe_su3_prod_su3(fixer[ivol],g,fixer[ivol]);
	    
	    if(rank==0 and eo==0 and ieo==0)
	      {
		su3 test;
		unsafe_su3_prod_su3_dag(test,fixer[ivol],fixer[ivol]);
		master_printf("test:\n");
		su3_print(test);
		master_printf("g %d:\n",ivol);
		su3_print(fixer[ivol]);
	      }

	    //lower external border must be sync.ed with upper internal border of lower node
	    //  upper internal border with same parity must be sent using buf_up[mu][par]
	    //    ""     ""      ""   ""   opp.    "     "  "  recv using buf_up[mu][!par]
	    //  lower external   ""   ""   same    "     "  "  recv using buf_dw[mu][!par]
	    //    ""     ""      ""   ""   opp.    "     "  "  sent using buf_dw[mu][par]
	    for(int mu=0;mu<NDIM;mu++)
	      {
		int f=loclx_neighup[ivol][mu];
		int b=loclx_neighdw[ivol][mu];
		if(f>=loc_vol) su3_copy(((su3*)send_buf)[loceo_of_loclx[f]-loc_volh],conf[ivol][mu]);
		if(b>=loc_vol) su3_copy(((su3*)send_buf)[loceo_of_loclx[b]-loc_volh],conf[b][mu]);
	      }
	  }
	THREAD_BARRIER();
	
	//communicate
	comm_start(eo_su3_comm);
	comm_wait(eo_su3_comm);
	
	//read out
	NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
	  if(loclx_parity[ivol]!=eo)
	    for(int mu=0;mu<NDIM;mu++)
	      {
		int f=loclx_neighup[ivol][mu];
		int b=loclx_neighdw[ivol][mu];
		if(f>=loc_vol) su3_copy(conf[ivol][mu],((su3*)recv_buf)[loceo_of_loclx[f]-loc_volh]);
		if(b>=loc_vol) su3_copy(conf[b][mu],((su3*)recv_buf)[loceo_of_loclx[b]-loc_volh]);
	      }
	THREAD_BARRIER();
      }
    
    set_borders_invalid(conf);
  }
  
  //put the Fourier Acceleration kernel of eq.3.6 of C.Davies paper
  void Fourier_accelerate_derivative(su3 *der)
  {
    GET_THREAD_ID();
    
    //su3_print(der[0]);
    
    //Fourier Transform
    fft4d(der,-1,1);
    
    //compute 4*\sum_mu sin^2(2*pi*(L_mu-1))
    double num=16;
    
    //put the kernel and the prefactor
    NISSA_PARALLEL_LOOP(imom,0,loc_vol)
      {
	//compute 4*\sum_mu sin^2(2*pi*ip_mu)
	double den=0;
	for(int mu=0;mu<NDIM;mu++)
	  {
	    double p=2*M_PI*glb_coord_of_loclx[imom][mu]/glb_size[mu];
	    den+=sqr(sin(0.5*p));
	  }
	den*=4;
	
	//master_printf("%lg %lg rt\n",den,num);
	// master_printf("%lg %lg re\n",den,su3_norm2(der[imom]));
	
	//build the factor
	double fact=num;
	if(fabs(den)>1e-10) fact/=den;
	
	//master_printf("fact %d %lg\n",imom,fact);
	
	//put the factor
	su3_prodassign_double(der[imom],fact);
      }
    THREAD_BARRIER();
    
    //Anti-Fourier Transform
    fft4d(der,+1,0);
    
    // su3_print(der[0]);
    
    //crash("");
  }
  
  //do all the fixing with Fourier acceleration
  void Landau_or_Coulomb_gauge_fix_FACC(quad_su3 *fixed_conf,su3 *fixer,int start_mu,double alpha,quad_su3 *ori_conf)
  {
    GET_THREAD_ID();
    
    su3 *der=nissa_malloc("der",loc_vol,su3);
    
    //take the derivative
    communicate_lx_quad_su3_borders(fixed_conf);
    NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
      {
	su3 temp;
	compute_Landau_or_Coulomb_functional_der(temp,fixed_conf,ivol,start_mu);
	unsafe_su3_traceless_anti_hermitian_part(der[ivol],temp);
	
	su3_prodassign_double(der[ivol],-0.5*alpha);
      }
    THREAD_BARRIER();
    
    //put the kernel
    Fourier_accelerate_derivative(der);
    
    //take the exponent
    NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
      {
	su3 temp;
	safe_anti_hermitian_exact_exponentiate(temp,der[ivol]);
	safe_su3_prod_su3(fixer[ivol],temp,fixer[ivol]);
      }
    set_borders_invalid(fixer);
    
    //transform
    gauge_transform_conf(fixed_conf,fixer,ori_conf);
    
    nissa_free(der);
  }
  
  //check if gauge fixed or not
  bool check_Landau_or_Coulomb_gauge_fixed(double &prec,double &func,quad_su3 *fixed_conf,int start_mu,double target_prec)
  {
    prec=compute_Landau_or_Coulomb_gauge_fixing_quality(fixed_conf,start_mu);
    func=compute_Landau_or_Coulomb_functional(fixed_conf,start_mu);
    bool get_out=(prec<=target_prec);
    return get_out;
  }
  
  THREADABLE_FUNCTION_4ARG(Landau_or_Coulomb_gauge_fix, quad_su3*,fixed_conf, quad_su3*,ext_conf, int,start_mu, double,target_prec)
  {
    GET_THREAD_ID();
    double time=-take_time();
    
    //store input conf if equal
    quad_su3 *ori_conf=ext_conf;
    if(fixed_conf==ori_conf)
      {
	ori_conf=nissa_malloc("ori_conf",loc_vol+bord_vol+edge_vol,quad_su3);
	vector_copy(ori_conf,ext_conf);
 	vector_copy(fixed_conf,ext_conf);
      }
    
    //fix overrelax probability
    //const double over_relax_prob=0.9;
    
    //fixing transformation
    su3 *fixer=nissa_malloc("fixer",loc_vol+bord_vol,su3);
    NISSA_LOC_VOL_LOOP(ivol)
      su3_put_to_id(fixer[ivol]);
    set_borders_invalid(fixer);
    
    int macro_iter=0,nmax_macro_iter=100;
    double prec,func,old_func;
    double alpha=0.08,retuning_fact=1.2;
    bool really_get_out=check_Landau_or_Coulomb_gauge_fixed(prec,func,fixed_conf,start_mu,target_prec);
    do
      {
	//go on fixing until reaching precision, or exceeding the
	//iteration count
        int iter=0,nmax_iter=100;
	bool get_out=false;
	do
	  {
	    master_printf("iter %d, quality: %16.16lg, functional: %16.16lg\n",iter,prec,func);
	    
	    //store
	    //old_prec=prec;
	    old_func=func;
	    
	    Landau_or_Coulomb_gauge_fix_FACC(fixed_conf,fixer,start_mu,alpha,ori_conf);
	    
	    //Landau_or_Coulomb_gauge_fix(fixed_conf,fixer,start_mu,over_relax_prob);
	    iter++;
	    
	    //print out the precision reached and the functional
	    get_out=check_Landau_or_Coulomb_gauge_fixed(prec,func,fixed_conf,start_mu,target_prec);
	    
	    //retune alpha
	    if(func>old_func) alpha/=sqr(retuning_fact);
	    else              alpha*=retuning_fact;
	    master_printf("Changing alpha to %lg\n",alpha);
	  }
	while(iter<nmax_iter and not get_out);
	
	//now we put the fixer on su3, and make a real transformation
	//on the basis of what we managed to fix
	THREAD_BARRIER();
	if(rank==0 && thread_id==0)
	  {
	    su3 test;
	    safe_su3_prod_su3_dag(test,fixer[0],fixer[0]);
	    master_printf("test:\n");
	    su3_print(test);
	    su3_summ_real(test,test,-1);
	    master_printf(" unitarity before: %.16lg\n",sqrt(su3_norm2(test)));
	    master_printf("g %d:\n",0);
	    su3_print(fixer[0]);
	  }
	THREAD_BARRIER();
	NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
	  su3_unitarize_maximal_trace_projecting(fixer[ivol]);
	set_borders_invalid(fixer);
	THREAD_BARRIER();
	if(rank==0 && thread_id==0)
	  {
	    su3 test;
	    safe_su3_prod_su3_dag(test,fixer[0],fixer[0]);
	    master_printf("test:\n");
	    su3_print(test);
	    su3_summ_real(test,test,-1);
	    master_printf(" unitarity after: %.16lg\n",sqrt(su3_norm2(test)));
	    master_printf("g %d:\n",0);
	    su3_print(fixer[0]);
	  }
	THREAD_BARRIER();
	gauge_transform_conf(fixed_conf,fixer,ori_conf);
	
	//check if really get out
	really_get_out=check_Landau_or_Coulomb_gauge_fixed(prec,func,fixed_conf,start_mu,target_prec);
	master_printf("macro-iter %d, quality: %16.16lg, functional: %16.16lg\n",macro_iter,prec,func);
	macro_iter++;
      }
    while(macro_iter<nmax_macro_iter and not really_get_out);
    
    //crash if this did not work
    if(not really_get_out) crash("unable to fix to precision %16.16lg in %d macro-iterations",target_prec,macro_iter);
    
    if(fixed_conf==ori_conf) nissa_free(ori_conf);
    
    master_printf("Gauge fix time: %lg\n",time+take_time());
  }
  THREADABLE_FUNCTION_END
  
  //wrappers
  void Landau_gauge_fix(quad_su3 *conf_out,quad_su3 *conf_in,double precision)
  {Landau_or_Coulomb_gauge_fix(conf_out,conf_in,0,precision);}
  void Coulomb_gauge_fix(quad_su3 *conf_out,quad_su3 *conf_in,double precision)
  {Landau_or_Coulomb_gauge_fix(conf_out,conf_in,1,precision);}
  
  //perform a random gauge transformation
  THREADABLE_FUNCTION_2ARG(perform_random_gauge_transform, quad_su3*,conf_out, quad_su3*,conf_in)
  {
    GET_THREAD_ID();
    
    //allocate fixing matrix
    su3 *fixm=nissa_malloc("fixm",loc_vol+bord_vol,su3);
    
    //extract random SU(3) matrix
    NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
      su3_put_to_rnd(fixm[ivol],loc_rnd_gen[ivol]);
    set_borders_invalid(fixm);
    
    //apply the transformation
    gauge_transform_conf(conf_out,fixm,conf_in);
    
    //free fixing matrix
    nissa_free(fixm);
  }
  THREADABLE_FUNCTION_END
  
  THREADABLE_FUNCTION_2ARG(perform_random_gauge_transform, quad_su3**,conf_out, quad_su3**,conf_in)
  {
    GET_THREAD_ID();
    
    //allocate fixing matrix
    su3 *fixm[2]={nissa_malloc("fixm_e",loc_volh+bord_volh,su3),nissa_malloc("fixm_o",loc_volh+bord_volh,su3)};
    
    //extract random SU(3) matrix
    NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
      su3_put_to_rnd(fixm[loclx_parity[ivol]][loceo_of_loclx[ivol]],loc_rnd_gen[ivol]);
    for(int eo=0;eo<2;eo++) set_borders_invalid(fixm[eo]);
    
    //apply the transformation
    gauge_transform_conf(conf_out,fixm,conf_in);
    
    //free fixing matrix
    for(int eo=0;eo<2;eo++) nissa_free(fixm[eo]);
  }
  THREADABLE_FUNCTION_END
}
