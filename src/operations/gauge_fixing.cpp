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

#include "gauge_fixing.hpp"

#define VERBOSITY_MASTER_PRINTF verbosity_lv1_master_printf
//#define VERBOSITY_MASTER_PRINTF verbosity_lv3_master_printf

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
  THREADABLE_FUNCTION_3ARG(gauge_transform_conf, quad_su3*,uout, su3*,g, const quad_su3*,uin)
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
      }
    THREAD_BARRIER();
    
    //global reduction
    double omega;
    double_vector_glb_collapse(&omega,loc_omega,loc_vol);
    nissa_free(loc_omega);
    
    return omega/glb_vol/NCOL;
  }
  
  //do all the fixing
  void Landau_or_Coulomb_gauge_fixing_overrelax(quad_su3 *fixed_conf,LC_gauge_fixing_pars_t::gauge_t gauge,double overrelax_prob,su3 *fixer,quad_su3 *ori_conf)
  {
    GET_THREAD_ID();
    
    for(int eo=0;eo<2;eo++)
      {
	NISSA_PARALLEL_LOOP(ieo,0,loc_volh)
	  {
	    int ivol=loclx_of_loceo[eo][ieo];
	    
	    //compute the derivative
	    su3 temp;
	    compute_Landau_or_Coulomb_functional_der(temp,fixed_conf,ivol,gauge);
	    
	    //dagger
	    su3 ref;
	    unsafe_su3_hermitian(ref,temp);
	    
	    //find the link that maximizes the trace
	    su3 g;
	    su3_unitarize_maximal_trace_projecting(g,ref);
	    
	    //square probabilistically
	    double p=rnd_get_unif(loc_rnd_gen+ivol,0,1);
	    if(p<overrelax_prob) safe_su3_prod_su3(g,g,g);
	    
	    //transform
	    local_gauge_transform(fixed_conf,g,ivol);
	    
	    //store the change
	    safe_su3_prod_su3(fixer[ivol],g,fixer[ivol]);
	    
	    //lower external border must be sync.ed with upper internal border of lower node
	    //  upper internal border with same parity must be sent using buf_up[mu][par]
	    //    ""     ""      ""   ""   opp.    "     "  "  recv using buf_up[mu][!par]
	    //  lower external   ""   ""   same    "     "  "  recv using buf_dw[mu][!par]
	    //    ""     ""      ""   ""   opp.    "     "  "  sent using buf_dw[mu][par]
	    for(int mu=0;mu<NDIM;mu++)
	      {
		int f=loclx_neighup[ivol][mu];
		int b=loclx_neighdw[ivol][mu];
		if(f>=loc_vol) su3_copy(((su3*)send_buf)[loceo_of_loclx[f]-loc_volh],fixed_conf[ivol][mu]);
		if(b>=loc_vol) su3_copy(((su3*)send_buf)[loceo_of_loclx[b]-loc_volh],fixed_conf[b][mu]);
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
		if(f>=loc_vol) su3_copy(fixed_conf[ivol][mu],((su3*)recv_buf)[loceo_of_loclx[f]-loc_volh]);
		if(b>=loc_vol) su3_copy(fixed_conf[b][mu],((su3*)recv_buf)[loceo_of_loclx[b]-loc_volh]);
	      }
	THREAD_BARRIER();
      }
    
    set_borders_invalid(fixed_conf);
  }
  
  //put the Fourier Acceleration kernel of eq.3.6 of C.Davies paper
  void Fourier_accelerate_derivative(su3 *der)
  {
    GET_THREAD_ID();
    
    //Fourier Transform
    fft4d(der,FFT_MINUS,FFT_NORMALIZE);
    
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
	
	//build the factor
	double fact=num;
	if(fabs(den)>1e-10) fact/=den;
	
	//put the factor
	su3_prodassign_double(der[imom],fact);
      }
    THREAD_BARRIER();
    
    //Anti-Fourier Transform
    fft4d(der,FFT_PLUS,FFT_NO_NORMALIZE);
  }
  
  //take exp(-0.5*alpha*der)
  void exp_der_alpha_half(su3 *g,su3 *der,double alpha)
  {
    GET_THREAD_ID();
    NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
      {
	su3 temp;
	su3_prod_double(temp,der[ivol],-0.5*alpha);
	safe_anti_hermitian_exact_exponentiate(g[ivol],temp);
      }
    set_borders_invalid(g);
  }
  
  //add the current fixer to previous one
  void add_current_transformation(su3 *fixer_out,const su3 *g,const su3 *fixer_in)
  {
    GET_THREAD_ID();
    
    //add current transformation
    NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
      safe_su3_prod_su3(fixer_out[ivol],g[ivol],fixer_in[ivol]);
    set_borders_invalid(fixer_out);
  }
  
  //adapt the value of alpha to minimize the functional
  void adapt_alpha(quad_su3 *fixed_conf,su3 *fixer,int start_mu,su3 *der,double &alpha,double alpha_def,quad_su3 *ori_conf,const double func_0,bool &use_adapt)
  {
#ifndef REPRODUCIBLE_RUN
    crash("need reproducible run to enable adaptative search");
#endif
    
    //store original fixer
    su3 *ori_fixer=nissa_malloc("ori_fixer",loc_vol+bord_vol,su3);
    vector_copy(ori_fixer,fixer);
    
    //current transform
    su3 *g=nissa_malloc("g",loc_vol,su3);
    
    int zero_curv,pos_curv,brack_vert;
    //int nneg_pos_vert=0;
    int iter=0;
    const int nadapt_iter_max=5;
    do
      {
	VERBOSITY_MASTER_PRINTF("---iter %d---\n",iter);
	//take the exponent
	exp_der_alpha_half(g,der,alpha);
	
	//compute three points at 0,alpha and 2alpha
	double F[3];
	F[0]=func_0;
	vector_copy(fixer,ori_fixer);
	// master_printf("Check: %lg %lg\n",func_0,compute_Landau_or_Coulomb_functional(fixed_conf,start_mu));
	
	for(int i=1;i<=2;i++)
	  {
	    add_current_transformation(fixer,g,fixer);
	    
	    //transform and compute potential
	    gauge_transform_conf(fixed_conf,fixer,ori_conf);
	    F[i]=compute_Landau_or_Coulomb_functional(fixed_conf,start_mu);
	  }
	
	//subtract 0
	for(int i=2;i>=0;i--) F[i]-=F[0];
	
	double c=F[0];
	double b=(4*F[1]-F[2]-3*F[0])/(2*alpha);
	double a=(F[2]-2*F[1]+F[0])/(2*sqr(alpha));
	
	VERBOSITY_MASTER_PRINTF("F:   %.16lg %.16lg %.16lg\n",F[0],F[1],F[2]);
	VERBOSITY_MASTER_PRINTF("abc: %lg %lg %lg\n",a,b,c);
	
	double vert=-b/(2*a);
	pos_curv=(a>0);
	zero_curv=(fabs(a)<1e-14);
	brack_vert=(2*alpha>vert);
	
	VERBOSITY_MASTER_PRINTF("Vertex position: %lg\n",vert);
	VERBOSITY_MASTER_PRINTF("Curvature is positive: %d\n",pos_curv);
	if(zero_curv)
	  {
	    VERBOSITY_MASTER_PRINTF("Curvature is compatible with zero (%lg), switching off temporarily the adaptative search\n",a);
	    alpha=alpha_def;
	  }
	VERBOSITY_MASTER_PRINTF("Bracketing the vertex: %d\n",brack_vert);
	if(not pos_curv)
	  {
	    alpha/=2.05897683269763;
	    VERBOSITY_MASTER_PRINTF("Decreasing alpha to %lg\n",alpha);
	  }
	else
	  {
	    if(not brack_vert) VERBOSITY_MASTER_PRINTF("Not bracketing the vertex, increasing alpha to %lg\n",alpha);
	    else               VERBOSITY_MASTER_PRINTF("Good, jumping to %lg\n",vert);
	  }
	alpha=vert;
	
	if(iter>=nadapt_iter_max)
	  {
	    VERBOSITY_MASTER_PRINTF("%d adaptative searches performed, switching temporarily off the adaptative search\n",iter);
	    alpha=alpha_def;
	  }
	
	// //compute average and stddev
	// double ave,dev;
	// ave_dev(ave,dev,F,3);
	// VERBOSITY_MASTER_PRINTF("F ave, dev: %lg %lg\n",ave,dev);
	// if(dev<fabs(ave)*1e-15)
	//   {
	//     master_printf("Switching off adaptative search\n");
	//     use_adapt=false;
	//   }
	
	iter++;
      }
    while(use_adapt and (iter<nadapt_iter_max) and not (zero_curv or (pos_curv and brack_vert)));
    
    //put back the fixer
    vector_copy(fixer,ori_fixer);
    nissa_free(ori_fixer);
    nissa_free(g);
  }
  
  //GCG stuff - taken from 1405.5812
  namespace GCG
  {
    su3 *prev_der;
    su3 *s;
    double *accum;
  }
  
  //allocate the stuff for GCG
  void allocate_GCG_stuff()
  {
    using namespace GCG;
    
    prev_der=nissa_malloc("prev_der",loc_vol,su3);
    s=nissa_malloc("s",loc_vol,su3);
    accum=nissa_malloc("accum",loc_vol,double);
    
    vector_reset(prev_der);
    vector_reset(s);
  }
  
  //free the stuff for GCG
  void free_GCG_stuff()
  {
    using namespace GCG;
    
    nissa_free(prev_der);
    nissa_free(s);
    nissa_free(accum);
  }
  
  //apply GCG acceleration
  void GCG_improve_gauge_fixer(su3 *der,bool &use_GCG,int iter)
  {
    using namespace GCG;
    
    GET_THREAD_ID();
    
    double beta;
    if(iter>1)
      {
	//denominator
	NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
	  accum[ivol]=real_part_of_trace_su3_prod_su3_dag(prev_der[ivol],prev_der[ivol]);
	THREAD_BARRIER();
	double den;
	double_vector_glb_collapse(&den,accum,loc_vol);
	VERBOSITY_MASTER_PRINTF("den: %lg\n",den);
	
	//numerator
	NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
	  if(den>1e-6) //means that |DA|<1e-12
	    {
	      su3 temp;
	      su3_subt(temp,der[ivol],prev_der[ivol]);
	      accum[ivol]=real_part_of_trace_su3_prod_su3_dag(der[ivol],temp);
	    }
	  else accum[ivol]=real_part_of_trace_su3_prod_su3_dag(der[ivol],der[ivol]);
	THREAD_BARRIER();
	double num;
	double_vector_glb_collapse(&num,accum,loc_vol);
	VERBOSITY_MASTER_PRINTF("num: %lg\n",num);
	
	//compute beta
	beta=num/den;
	if(beta<0) beta=0;
	
	//switch off beta if even smaller
	const double gcg_tol=1e-20;
	if(fabs(num)<gcg_tol or fabs(den)<gcg_tol)
	  {
	    beta=0;
	    use_GCG=false;
	    VERBOSITY_MASTER_PRINTF("Switching off GCG at iter %d, fabs(num)[%lg]<gcg_tol[%lg] or fabs(den)[%lg]<gcg_tol[%lg]\n",iter,fabs(num),gcg_tol,fabs(den),gcg_tol);
	  }
      }
    else beta=0;
    VERBOSITY_MASTER_PRINTF("beta: %lg\n",beta);
    
    //store prev_der, increase s (der) and store prev_s
    vector_copy(prev_der,der);
    if(iter==0) vector_copy(s,der);
    else        double_vector_summ_double_vector_prod_double((double*)s,(double*)der,(double*)s,beta,loc_vol*sizeof(su3)/sizeof(double));
  }
  
  //do all the fixing exponentiating
  void Landau_or_Coulomb_gauge_fixing_exponentiate(quad_su3 *fixed_conf,su3 *fixer,LC_gauge_fixing_pars_t::gauge_t gauge,
						   double &alpha,double alpha_def,quad_su3 *ori_conf,const double func_0,const bool &use_FACC,bool &use_adapt,bool &use_GCG,int iter)
  {
    using namespace GCG;
    
    GET_THREAD_ID();
    
    //take the derivative
    su3 *der=nissa_malloc("der",loc_vol,su3);
    communicate_lx_quad_su3_borders(fixed_conf);
    NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
      {
	su3 temp;
	compute_Landau_or_Coulomb_functional_der(temp,fixed_conf,ivol,gauge);
	unsafe_su3_traceless_anti_hermitian_part(der[ivol],temp);
      }
    THREAD_BARRIER();
    
    //put the kernel
    if(use_FACC) Fourier_accelerate_derivative(der);
    
    //make the CG improvement
    if(use_GCG) GCG_improve_gauge_fixer(der,use_GCG,iter);
    
    //decides what to use
    su3 *v=(use_GCG?s:der);
    
    //take the exponent with alpha
    su3 *g=nissa_malloc("g",loc_vol,su3);
    
    if(use_adapt) adapt_alpha(fixed_conf,fixer,gauge,v,alpha,alpha_def,ori_conf,func_0,use_adapt);
    exp_der_alpha_half(g,v,alpha);
    
    //put the transformation
    add_current_transformation(fixer,g,fixer);
    gauge_transform_conf(fixed_conf,fixer,ori_conf);
    
    nissa_free(der);
    nissa_free(g);
  }
  
  //check if gauge fixed or not
  bool check_Landau_or_Coulomb_gauge_fixed(double &prec,double &func,quad_su3 *fixed_conf,LC_gauge_fixing_pars_t::gauge_t gauge,double target_prec)
  {
    prec=compute_Landau_or_Coulomb_gauge_fixing_quality(fixed_conf,gauge);
    func=compute_Landau_or_Coulomb_functional(fixed_conf,gauge);
    bool get_out=(prec<=target_prec);
    return get_out;
  }
  
  THREADABLE_FUNCTION_3ARG(Landau_or_Coulomb_gauge_fix, quad_su3*,fixed_conf, LC_gauge_fixing_pars_t*,pars, quad_su3*,ext_conf)
  {
    GET_THREAD_ID();
    double time=-take_time();
    
    const bool use_fft_acc=pars->use_fft_acc;
    bool use_adapt=pars->use_adaptative_search;
    bool use_GCG=pars->use_generalized_cg;
    if(pars->use_generalized_cg) allocate_GCG_stuff();
    
    //store input conf if equal
    quad_su3 *ori_conf=ext_conf;
    if(fixed_conf==ori_conf)
      {
	ori_conf=nissa_malloc("ori_conf",loc_vol+bord_vol+edge_vol,quad_su3);
	vector_copy(ori_conf,ext_conf);
      }
    vector_copy(fixed_conf,ext_conf);
    
    //fixing transformation
    su3 *fixer=nissa_malloc("fixer",loc_vol+bord_vol,su3);
    NISSA_LOC_VOL_LOOP(ivol)
      su3_put_to_id(fixer[ivol]);
    set_borders_invalid(fixer);
    
    double prec,func,old_func;
    double alpha=pars->alpha_exp;
    bool really_get_out=check_Landau_or_Coulomb_gauge_fixed(prec,func,fixed_conf,pars->gauge,pars->target_precision);
    int iter=0;
    do
      {
	//go on fixing until reaching precision, or exceeding the iteration count
	bool get_out=false;
	do
	  {
	    master_printf("iter: %d quality: %16.16lg functional: %16.16lg\n",iter,prec,func);
	    
	    //store
	    old_func=func;
	    
	    switch(pars->method)
	      {
	      case LC_gauge_fixing_pars_t::exponentiate:
		Landau_or_Coulomb_gauge_fixing_exponentiate(fixed_conf,fixer,pars->gauge,alpha,pars->alpha_exp,ori_conf,old_func,use_fft_acc,use_adapt,use_GCG,iter);break;
	      case LC_gauge_fixing_pars_t::overrelax:
		Landau_or_Coulomb_gauge_fixing_overrelax(fixed_conf,pars->gauge,pars->overrelax_prob,fixer,ori_conf);break;
	      default:
		crash("unknown method %d",pars->method);
	      }
	    iter++;
	    
	    //print out the precision reached and the functional
	    get_out=false;
	    get_out|=check_Landau_or_Coulomb_gauge_fixed(prec,func,fixed_conf,pars->gauge,pars->target_precision);
	    get_out|=(not (iter<pars->nmax_iterations));
	    get_out|=(not (iter%pars->unitarize_each==0));
	    
	    //switch off adaptative search if precision is too small
	    const double adapt_tol=1e-13;
	    if(use_adapt and prec<adapt_tol)
	      {
		master_printf("Switching off adaptative search\n");
		use_adapt=false;
	      }
	  }
	while(not get_out);
	
	//now we put the fixer on su3, and make a real transformation
	//on the basis of what we managed to fix
	//THREAD_BARRIER();
	// if(rank==0 && thread_id==0)
	//   {
	//     su3 test;
	//     safe_su3_prod_su3_dag(test,fixer[0],fixer[0]);
	//     master_printf("test:\n");
	//     su3_print(test);
	//     su3_summ_real(test,test,-1);
	//     master_printf(" unitarity before: %.16lg\n",sqrt(su3_norm2(test)));
	//     master_printf("g %d:\n",0);
	//     su3_print(fixer[0]);
	//   }
	// THREAD_BARRIER();
	NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
	  su3_unitarize_explicitly_inverting(fixer[ivol],fixer[ivol]);
	set_borders_invalid(fixer);
	// THREAD_BARRIER();
	// if(rank==0 && thread_id==0)
	//   {
	//     su3 test;
	//     safe_su3_prod_su3_dag(test,fixer[0],fixer[0]);
	//     master_printf("test:\n");
	//     su3_print(test);
	//     su3_summ_real(test,test,-1);
	//     master_printf(" unitarity after: %.16lg\n",sqrt(su3_norm2(test)));
	//     master_printf("g %d:\n",0);
	//     su3_print(fixer[0]);
	//   }
	// THREAD_BARRIER();
	gauge_transform_conf(fixed_conf,fixer,ori_conf);
	
	//check if really get out
	really_get_out=false;
	really_get_out|=check_Landau_or_Coulomb_gauge_fixed(prec,func,fixed_conf,pars->gauge,pars->target_precision);
	really_get_out|=(not (iter<pars->nmax_iterations));
      }
    while(not really_get_out);
    
    //crash if this did not work
    if(not really_get_out) crash("unable to fix to precision %16.16lg in %d iterations",pars->target_precision,pars->nmax_iterations);
    
    //free
    if(ori_conf!=ext_conf) nissa_free(ori_conf);
    nissa_free(fixer);
    if(pars->use_generalized_cg) free_GCG_stuff();
    
    master_printf("Gauge fix time: %lg\n",time+take_time());
  }
  THREADABLE_FUNCTION_END
  
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
