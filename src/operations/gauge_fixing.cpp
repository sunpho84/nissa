#ifdef HAVE_CONFIG_H
# include "config.hpp"
#endif

#include <math.h>
#include <string.h>

#include <base/debug.hpp>
#include <base/field.hpp>
#include <base/random.hpp>
#include <base/vectors.hpp>
#include <communicate/borders.hpp>
#include <geometry/geometry_eo.hpp>
#include <geometry/geometry_lx.hpp>
#include <linalgs/linalgs.hpp>
#include <linalgs/reduce.hpp>
#include <new_types/complex.hpp>
#include <new_types/su3_op.hpp>
#include <operations/fft.hpp>
#include <operations/gauge_fixing.hpp>
#include <routines/ios.hpp>
#include <routines/mpi_routines.hpp>
#include <threads/threads.hpp>

#define VERBOSITY_MASTER_PRINTF verbosity_lv1_master_printf
//#define VERBOSITY_MASTER_PRINTF verbosity_lv3_master_printf

namespace nissa
{
  //apply the passed transformation to the point
  CUDA_HOST_AND_DEVICE void local_gauge_transform(quad_su3 *conf,su3 g,int ivol)
  {
    // for each dir...
    for(int mu=0;mu<NDIM;mu++)
      {
        int b=loclxNeighdw[ivol][mu];
        
        //perform local gauge transform
        safe_su3_prod_su3(conf[ivol][mu],g,conf[ivol][mu]);
        safe_su3_prod_su3_dag(conf[b][mu],conf[b][mu],g);
      }
  }
  
  void gauge_transform_conf(LxField<quad_su3>& uout,
			    const LxField<su3>& g,
			    const LxField<quad_su3>& uin)
  {
    g.updateHalo();
    
    //transform
    PAR(0,locVol,
	CAPTURE(TO_WRITE(uout),
		TO_READ(uin),
		TO_READ(g)),
	ivol,
	{
	  for(int mu=0;mu<NDIM;mu++)
	    {
	      su3 temp;
	      unsafe_su3_prod_su3_dag(temp,uin[ivol][mu],g[loclxNeighup[ivol][mu]]);
	      unsafe_su3_prod_su3(uout[ivol][mu],g[ivol],temp);
	    }
	});
  }
  
  //e/o version
  void gauge_transform_conf(EoField<quad_su3>& uout,
			    const EoField<su3>& g,
			    const EoField<quad_su3>& uin)
  {
    //communicate borders
    g.updateHalo();
    
    //transform
    for(int par=0;par<2;par++)
      PAR(0,locVolh,
	  CAPTURE(par,
		  TO_WRITE(uout),
		  TO_READ(g),
		  TO_READ(uin)),
	  ivol,
	  {
	    for(int mu=0;mu<NDIM;mu++)
	      {
		su3 temp;
		unsafe_su3_prod_su3_dag(temp,uin[par][ivol][mu],g[!par][loceo_neighup[par][ivol][mu]]);
		unsafe_su3_prod_su3(uout[par][ivol][mu],g[par][ivol],temp);
	      }
	  });
  }
  
  //transform a color field
  void gauge_transform_color(eo_ptr<color> out,eo_ptr<su3> g,eo_ptr<color> in)
  {
    crash("reimplement");
    // //communicate borders
    // communicate_ev_and_od_su3_borders(g);
    
    // //transform
    // for(int par=0;par<2;par++)
    //   NISSA_PARALLEL_LOOP(ivol,0,locVolh)
    // 	safe_su3_prod_color(out[par][ivol],g[par][ivol],in[par][ivol]);
    // NISSA_PARALLEL_LOOP_END;
    
    // //invalidate borders
    // set_borders_invalid(out[0]);
    // set_borders_invalid(out[1]);
  }
  
  //determine the gauge transformation bringing to temporal gauge with T-1 timeslice diferent from id
  void find_temporal_gauge_fixing_matr(su3 *fixm,quad_su3 *u)
  {
    crash("reimplement");
    // int loc_slice_area=locSize[1]*locSize[2]*locSize[3];
    // su3 *buf=NULL;
    
    // //if the number of ranks in the 0 dir is greater than 1 allocate room for border
    // if(nrank_dir[0]>1) buf=nissa_malloc("buf",loc_slice_area,su3);
    
    // //if we are on first rank slice put to identity the t=0 slice, otherwise receive it from previous rank slice
    // if(rank_coord[0]==0)
    //   {
    // 	NISSA_LOC_VOL_LOOP(ivol)
    // 	  if(glbCoordOfLoclx[ivol][0]==0)
    // 	    su3_put_to_id(fixm[ivol]);
    //   }
    // else
    //   if(nrank_dir[0]>1)
    // 	MPI_Recv((void*)fixm,loc_slice_area,MPI_SU3,rank_neighdw[0],252,cart_comm,MPI_STATUS_IGNORE);
    
    // //now go ahead along t
    // Coords c;
    // //loop over spatial slice
    // for(c[1]=0;c[1]<locSize[1];c[1]++)
    //   for(c[2]=0;c[2]<locSize[2];c[2]++)
    // 	for(c[3]=0;c[3]<locSize[3];c[3]++)
    // 	  {
    // 	    //bulk
    // 	    for(c[0]=1;c[0]<locSize[0];c[0]++)
    // 	      {
    // 		int icurr=loclx_of_coord(c);
    // 		c[0]--;int iback=loclx_of_coord(c);c[0]++;
		
    // 		unsafe_su3_prod_su3(fixm[icurr],fixm[iback],u[iback][0]);
    // 	      }
    // 	    //border
    // 	    if(nrank_dir[0]>1)
    // 	      {
    // 		c[0]=locSize[0]-1;int iback=loclx_of_coord(c);
    // 		c[0]=0;int icurr=loclx_of_coord(c);
		
    // 		unsafe_su3_prod_su3(buf[icurr],fixm[iback],u[iback][0]);
    // 	      }
	    
    // 	  }
    
    // //if we are not on last slice of rank send g to next slice
    // if(rank_coord[0]!=(nrank_dir[0]-1) && nrank_dir[0]>1)
    //   MPI_Send((void*)buf,loc_slice_area,MPI_SU3,rank_neighup[0],252,cart_comm);
    
    // if(nrank_dir[0]>1) nissa_free(buf);
  }
  
  ////////////////////////////////////// Landau or Coulomb gauges ///////////////////////////////////////////////////////
  
  //compute the functional on a single point
  double compute_Landau_or_Coulomb_functional(const LxField<quad_su3>& conf,
					      const int& ivol,
					      const int& start_mu)
  {
    double F=0;
    
    for(int mu=start_mu;mu<NDIM;mu++)
      {
	F-=su3_real_trace(conf[ivol][mu]);
	F-=su3_real_trace(conf[loclxNeighdw[ivol][mu]][mu]);
      }
    
    return F;
  }
  
  //derivative of the functional
  template <typename C>
  CUDA_HOST_AND_DEVICE INLINE_FUNCTION
  void compute_Landau_or_Coulomb_functional_der(su3& out,
						const C& conf,
						const int& ivol,
						const int& start_mu)
  {
    su3_put_to_zero(out);
    
    for(int mu=start_mu;mu<NDIM;mu++)
      {
	su3_summassign(out,conf[ivol][mu]);
	su3_summassign_su3_dag(out,conf[loclxNeighdw[ivol][mu]][mu]);
      }
  }
  
  //compute the functional that gets minimised
  double compute_Landau_or_Coulomb_functional(const LxField<quad_su3>& conf,
					      const int& start_mu,
					      const LxField<double> *F_offset=nullptr,
					      LxField<double> *ext_loc_F=nullptr)
  {
    
    LxField<double> *_loc_F=ext_loc_F;
    if(ext_loc_F==nullptr)
      _loc_F=new LxField<double>("loc_F");
    
    LxField<double>& loc_F=*_loc_F;
    
    PAR(0,locVol,
	CAPTURE(F_offset,start_mu,
		TO_WRITE(loc_F),
		TO_READ(conf)),
	ivol,
	{
	  if(F_offset) loc_F[ivol]=-(*F_offset)[ivol];
	  else         loc_F[ivol]=0;
	  
	  for(int mu=start_mu;mu<NDIM;mu++)
	    loc_F[ivol]-=su3_real_trace(conf[ivol][mu]);
	});
    
    //collapse
    double F;
    glb_reduce(&F,loc_F,locVol);
    
    if(ext_loc_F==nullptr)
      delete _loc_F;
    
    return F;
  }
  
  //compute the quality of the Landau or Coulomb gauge fixing
  double compute_Landau_or_Coulomb_gauge_fixing_quality(const LxField<quad_su3>& conf,
							const int& start_mu)
  {
    conf.updateHalo();
    
    LxField<double> loc_omega("loc_omega");
    
    PAR(0,locVol,
	CAPTURE(start_mu,
		TO_WRITE(loc_omega),
		TO_READ(conf)),
	ivol,
	{
	  su3 delta;
	  su3_put_to_zero(delta);
	  
	  for(int mu=start_mu;mu<NDIM;mu++)
	    {
	      su3_subtassign(delta,conf[ivol][mu]);
	      su3_summassign(delta,conf[loclxNeighdw[ivol][mu]][mu]);
	    }
	  
	  //take 2 the traceless anti-hermitian part
	  su3 delta_TA;
	  unsafe_su3_traceless_anti_hermitian_part(delta_TA,delta);
	  loc_omega[ivol]=4*su3_norm2(delta_TA);
	});
    
    //global reduction
    double omega;
    glb_reduce(&omega,loc_omega,locVol);
    
    return omega/glbVol/NCOL;
  }
  
  //do all the fixing
  void Landau_or_Coulomb_gauge_fixing_overrelax(LxField<quad_su3>& fixed_conf,
						const LC_gauge_fixing_pars_t::gauge_t& gauge,
						const double& overrelax_prob,
						const LxField<su3>& fixer,
						const LxField<quad_su3>& ori_conf)
  {
    crash("reimplement");
    
    // for(int eo=0;eo<2;eo++)
    //   {
    // 	NISSA_PARALLEL_LOOP(ieo,0,locVolh)
    // 	  {
    // 	    int ivol=loclx_of_loceo[eo][ieo];
	    
    // 	    //compute the derivative
    // 	    su3 temp;
    // 	    compute_Landau_or_Coulomb_functional_der(temp,fixed_conf,ivol,gauge);
	    
    // 	    //dagger
    // 	    su3 ref;
    // 	    unsafe_su3_hermitian(ref,temp);
	    
    // 	    //find the link that maximizes the trace
    // 	    su3 g;
    // 	    su3_unitarize_maximal_trace_projecting(g,ref);
	    
    // 	    //square probabilistically
    // 	    double p=rnd_get_unif(loc_rnd_gen+ivol,0,1);
    // 	    if(p<overrelax_prob) safe_su3_prod_su3(g,g,g);
	    
    // 	    //transform
    // 	    local_gauge_transform(fixed_conf,g,ivol);
	    
    // 	    //store the change
    // 	    safe_su3_prod_su3(fixer[ivol],g,fixer[ivol]);
	    
    // 	    //lower external border must be sync.ed with upper internal border of lower node
    // 	    //  upper internal border with same parity must be sent using buf_up[mu][par]
    // 	    //    ""     ""      ""   ""   opp.    "     "  "  recv using buf_up[mu][!par]
    // 	    //  lower external   ""   ""   same    "     "  "  recv using buf_dw[mu][!par]
    // 	    //    ""     ""      ""   ""   opp.    "     "  "  sent using buf_dw[mu][par]
    // 	    for(int mu=0;mu<NDIM;mu++)
    // 	      {
    // 		int f=loclxNeighup[ivol][mu];
    // 		int b=loclxNeighdw[ivol][mu];
    // 		if(f>=locVol) su3_copy(((su3*)send_buf)[loceo_of_loclx[f]-locVolh],fixed_conf[ivol][mu]);
    // 		if(b>=locVol) su3_copy(((su3*)send_buf)[loceo_of_loclx[b]-locVolh],fixed_conf[b][mu]);
    // 	      }
    // 	  }
    // 	NISSA_PARALLEL_LOOP_END;
    // 	THREAD_BARRIER();
	
    // 	//communicate
    // 	comm_start(eo_su3_comm);
    // 	comm_wait(eo_su3_comm);
	
    // 	//read out
    // 	NISSA_PARALLEL_LOOP(ivol,0,locVol)
    // 	  if(loclx_parity[ivol]!=eo)
    // 	    for(int mu=0;mu<NDIM;mu++)
    // 	      {
    // 		int f=loclxNeighup[ivol][mu];
    // 		int b=loclxNeighdw[ivol][mu];
    // 		if(f>=locVol) su3_copy(fixed_conf[ivol][mu],((su3*)recv_buf)[loceo_of_loclx[f]-locVolh]);
    // 		if(b>=locVol) su3_copy(fixed_conf[b][mu],((su3*)recv_buf)[loceo_of_loclx[b]-locVolh]);
    // 	      }
    // 	NISSA_PARALLEL_LOOP_END;
    // 	THREAD_BARRIER();
    //   }
    
    // set_borders_invalid(fixed_conf);
  }
  
  //put the Fourier Acceleration kernel of eq.3.6 of C.Davies paper
  void Fourier_accelerate_derivative(LxField<su3>& der)
  {
    crash("reimplement");
    // //Fourier Transform
    // fft4d(der,FFT_MINUS,FFT_NORMALIZE);
    
    // //compute 4*\sum_mu sin^2(2*pi*(L_mu-1))
    // double num=16;
    
    // //put the kernel and the prefactor
    // NISSA_PARALLEL_LOOP(imom,0,locVol)
    //   {
    // 	//compute 4*\sum_mu sin^2(2*pi*ip_mu)
    // 	double den=0;
    // 	for(int mu=0;mu<NDIM;mu++)
    // 	  {
    // 	    double p=2*M_PI*glbCoordOfLoclx[imom][mu]/glbSize[mu];
    // 	    den+=sqr(sin(0.5*p));
    // 	  }
    // 	den*=4;
	
    // 	//build the factor
    // 	double fact=num;
    // 	if(fabs(den)>1e-10) fact/=den;
	
    // 	//put the factor
    // 	su3_prodassign_double(der[imom],fact);
    //   }
    // NISSA_PARALLEL_LOOP_END;
    // THREAD_BARRIER();
    
    // //Anti-Fourier Transform
    // fft4d(der,FFT_PLUS,FFT_NO_NORMALIZE);
  }
  
  //take exp(-0.5*alpha*der)
  void exp_der_alpha_half(LxField<su3>& g,
			  const LxField<su3>& der,
			  const double& alpha)
  {
    PAR(0,locVol,
	CAPTURE(alpha,
		TO_WRITE(g),
		TO_READ(der)),
	ivol,
	{
	  su3 temp;
	  su3_prod_double(temp,der[ivol],-0.5*alpha);
	  safe_anti_hermitian_exact_exponentiate(g[ivol],temp);
	});
  }
  
  /// Add the current fixer to previous one
  void add_current_transformation(LxField<su3>& fixer_out,
				  const LxField<su3>& g,
				  const LxField<su3>& fixer_in)
  {
    PAR(0,locVol,
	CAPTURE(TO_WRITE(fixer_out),
		TO_READ(fixer_in),
		TO_READ(g)),
	ivol,
	{
	  safe_su3_prod_su3(fixer_out[ivol],g[ivol],fixer_in[ivol]);
	});
  }
  
  //adapt the value of alpha to minimize the functional
  double adapt_alpha(LxField<quad_su3>& fixed_conf,
		     LxField<su3>& fixer,
		     const int& start_mu,
		     const LxField<su3>& der,
		     const double& alpha_def,
		     const LxField<quad_su3>& ori_conf,
		     const LxField<double> *F_offset,
		     const double& func,
		     const bool& use_adapt,
		     int& nskipped_adapt)
  {
#ifndef REPRODUCIBLE_RUN
    crash("need reproducible run to enable adaptative search");
#endif
    
    //first guess
    double alpha=alpha_def;
    
    //store original fixer
    LxField<su3> ori_fixer("ori_fixer",WITH_HALO);
    ori_fixer=fixer;
    
    //current transform
    LxField<su3> g("g");
    
    //int nneg_pos_vert=0;
    int iter=0;
    const int nadapt_iter_max=5;
    bool found=false,give_up=false;
    do
      {
	VERBOSITY_MASTER_PRINTF("---iter %d---\n",iter);
	//take the exponent
	exp_der_alpha_half(g,der,alpha);
	
	//compute three points at 0,alpha and 2alpha
	double F[3];
	F[0]=0;
	fixer=ori_fixer;
	// master_printf("Check: %lg %lg\n",func_0,compute_Landau_or_Coulomb_functional(fixed_conf,start_mu));
	
	for(int i=1;i<=2;i++)
	  {
	    add_current_transformation(fixer,g,fixer);
	    
	    //transform and compute potential
	    gauge_transform_conf(fixed_conf,fixer,ori_conf);
	    F[i]=compute_Landau_or_Coulomb_functional(fixed_conf,start_mu,F_offset);
	  }
	
	double c=F[0];
	double b=(4*F[1]-F[2]-3*F[0])/(2*alpha);
	double a=(F[2]-2*F[1]+F[0])/(2*sqr(alpha));
	
	VERBOSITY_MASTER_PRINTF("x:   %.16lg %.16lg %.16lg\n",0.0,alpha,2*alpha);
	VERBOSITY_MASTER_PRINTF("F:   %.16lg %.16lg %.16lg\n",F[0],F[1],F[2]);
	VERBOSITY_MASTER_PRINTF("abc: %lg %lg %lg\n",a,b,c);
	
	const double a_tol=1e-14;
	bool pos_curv=(a>a_tol);
	if(not pos_curv)
	  {
	    VERBOSITY_MASTER_PRINTF("Curvature %lg is not positive within tolerance (%lg), switching temporarily off the adaptative search\n",a,a_tol);
	    give_up=true;
	  }
	else
	  {
	    double vert=-b/(2*a);
	    pos_curv=(a>0);
	    bool brack_vert=(fabs(2*alpha)>fabs(vert));
	    
	    VERBOSITY_MASTER_PRINTF("Vertex position: %lg\n",vert);
	    VERBOSITY_MASTER_PRINTF("Curvature is positive: %d\n",pos_curv);
	    
	    alpha=vert;
	    if(not brack_vert) VERBOSITY_MASTER_PRINTF("Not bracketing the vertex, changing alpha to %lg\n",alpha);
	    else
	      {
		found=true;
		VERBOSITY_MASTER_PRINTF("Bracketting the vertex, jumping to %lg\n",alpha);
	      }
	  }
	
	iter++;
	
	//check that not too many iterations have been performed
	if(iter>=nadapt_iter_max)
	  {
	    VERBOSITY_MASTER_PRINTF("%d adaptative searches performed, switching off temporarily the adaptative search\n",iter);
	    give_up=true;
	  }
      }
    while((not found) and (not give_up));
    
    if(found) nskipped_adapt=0;
    else
      {
	nskipped_adapt++;
	alpha=alpha_def;
      }
    
    //put back the fixer
    *fixer=*ori_fixer;
    
    return alpha;
  }
  
  //GCG stuff - taken from 1405.5812
  namespace GCG
  {
    CUDA_MANAGED LxField<su3> *_prev_der;
    CUDA_MANAGED LxField<su3> *_s;
    CUDA_MANAGED LxField<double> *_accum;
  }
  
  //allocate the stuff for GCG
  void allocate_GCG_stuff()
  {
    using namespace GCG;
    
    _prev_der=new LxField<su3>("prev_der");
    _s=new LxField<su3>("s");
    _accum=new LxField<double>("accum");
    
    _prev_der->reset();
    _s->reset();
  }
  
  //free the stuff for GCG
  void free_GCG_stuff()
  {
    using namespace GCG;
    
    delete _prev_der;
    delete _s;
    delete _accum;
  }
  
  //apply GCG acceleration
  void GCG_improve_gauge_fixer(LxField<su3>& der,
			       bool &use_GCG,
			       const int& iter)
  {
    auto& prev_der=*GCG::_prev_der;
    auto& accum=*GCG::_accum;
    auto& s=*GCG::_s;
    
    double beta;
    if(iter>1)
      {
	//denominator
	PAR(0,locVol,
	    CAPTURE(TO_WRITE(accum),
		    TO_READ(prev_der)),
	    ivol,
	    {
	      accum[ivol]=real_part_of_trace_su3_prod_su3_dag(prev_der[ivol],prev_der[ivol]);
	    });
	
	double den;
	glb_reduce(&den,accum,locVol);
	VERBOSITY_MASTER_PRINTF("den: %lg\n",den);
	
	//numerator
	PAR(0,locVol,
	    CAPTURE(den,
		    TO_WRITE(accum),
		    TO_READ(prev_der),
		    TO_READ(der)),
	    ivol,
	    {
	      if(den>1e-6) //means that |DA|<1e-12
		{
		  su3 temp;
		  su3_subt(temp,der[ivol],prev_der[ivol]);
		  accum[ivol]=real_part_of_trace_su3_prod_su3_dag(der[ivol],temp);
		}
	      else accum[ivol]=real_part_of_trace_su3_prod_su3_dag(der[ivol],der[ivol]);
	    });
	
	double num;
	glb_reduce(&num,accum,locVol);
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
    prev_der=der;
    if(iter%10==0) s=der;
    else
      FOR_EACH_SITE_DEG_OF_FIELD(s,
				 CAPTURE(beta,
					 TO_READ(der),
					 TO_WRITE(s)),
				 ivol,iDeg,
      {
	auto& si=s(ivol,iDeg);
	si=der(ivol,iDeg)+si*beta;
      });
  }
  
  //do all the fixing exponentiating
  void Landau_or_Coulomb_gauge_fixing_exponentiate(LxField<quad_su3>& fixed_conf,
						   LxField<su3>& fixer,
						   const LC_gauge_fixing_pars_t::gauge_t& gauge,
						   const double alpha_def,
						   const LxField<quad_su3>& ori_conf,
						   const LxField<double>& F_offset,
						   const double& func,
						   const bool &use_FACC,
						   const bool &use_adapt,
						   int &nskipped_adapt,
						   bool &use_GCG,
						   const int& iter)
  {
    auto& s=*GCG::_s;
    
    fixed_conf.updateHalo();
    
    //take the derivative
    LxField<su3> der("der");
    
    PAR(0,locVol,
	CAPTURE(gauge,TO_WRITE(der),
		TO_READ(fixed_conf)),
	ivol,
      {
	su3 temp;
	compute_Landau_or_Coulomb_functional_der(temp,fixed_conf,ivol,gauge);
	unsafe_su3_traceless_anti_hermitian_part(der[ivol],temp);
      });
    
    //put the kernel
    if(use_FACC) Fourier_accelerate_derivative(der);
    
    //make the CG improvement
    if(use_GCG) GCG_improve_gauge_fixer(der,use_GCG,iter);
    
    //decides what to use
    LxField<su3>& v=(use_GCG?s:der);
    
    //take the exponent with alpha
    LxField<su3> g("g");
    
    //set alpha
    double alpha;
    if(use_adapt) alpha=adapt_alpha(fixed_conf,fixer,gauge,v,alpha_def,ori_conf,&F_offset,func,use_adapt,nskipped_adapt);
    else          alpha=alpha_def;
    
    exp_der_alpha_half(g,v,alpha);
    
    //put the transformation
    add_current_transformation(fixer,g,fixer);
    gauge_transform_conf(fixed_conf,fixer,ori_conf);
  }
  
  //check if gauge fixed or not
  bool check_Landau_or_Coulomb_gauge_fixed(double &prec,
					   double &func,
					   const LxField<quad_su3>& fixed_conf,
					   const LC_gauge_fixing_pars_t::gauge_t& gauge,
					   const double& target_prec,
					   LxField<double>* loc_F)
  {
    prec=compute_Landau_or_Coulomb_gauge_fixing_quality(fixed_conf,gauge);
    func=compute_Landau_or_Coulomb_functional(fixed_conf,gauge,nullptr,loc_F);
    
    const bool get_out=
      (prec<=target_prec);
    
    return get_out;
  }
  
  void Landau_or_Coulomb_gauge_fix(LxField<quad_su3>& fixed_conf,
				   const LC_gauge_fixing_pars_t& pars,
				   const LxField<quad_su3>& ext_conf)
  {
    if(fixed_conf==ext_conf)
      {
	LxField<quad_su3> tmp("tmp");
	tmp=ext_conf;
	Landau_or_Coulomb_gauge_fix(fixed_conf,pars,tmp);
      }
    else
      {
	double time=-take_time();
	
	const bool use_fft_acc=pars.use_fft_acc;
	bool use_adapt=pars.use_adaptative_search;
	bool use_GCG=pars.use_generalized_cg;
	if(pars.use_generalized_cg) allocate_GCG_stuff();
	
	fixed_conf=ext_conf;
	
	//fixing transformation
	LxField<su3> fixer("fixer",WITH_HALO);
	PAR(0,locVol,
	    CAPTURE(TO_WRITE(fixer)),
	    ivol,
	    {
	      su3_put_to_id(fixer[ivol]);
	    });
	
	LxField<double> F_offset("F_offset");
	double prec,func;
	bool really_get_out=check_Landau_or_Coulomb_gauge_fixed(prec,func,fixed_conf,pars.gauge,pars.target_precision,&F_offset);
	int iter=0,nskipped_adapt=0;
	do
	  {
	    //go on fixing until reaching precision, or exceeding the iteration count
	    bool get_out=false;
	    do
	      {
		master_printf("iter: %d quality: %16.16lg functional: %16.16lg\n",iter,prec,func);
		
		switch(pars.method)
		  {
		  case LC_gauge_fixing_pars_t::exponentiate:
		    Landau_or_Coulomb_gauge_fixing_exponentiate(fixed_conf,fixer,pars.gauge,pars.alpha_exp,ext_conf,F_offset,func,use_fft_acc,use_adapt,nskipped_adapt,use_GCG,iter);break;
		  case LC_gauge_fixing_pars_t::overrelax:
		    Landau_or_Coulomb_gauge_fixing_overrelax(fixed_conf,pars.gauge,pars.overrelax_prob,fixer,ext_conf);break;
		  default:
		    crash("unknown method %d",pars.method);
		  }
		iter++;
		
		//print out the precision reached and the functional
		get_out=false;
		get_out|=check_Landau_or_Coulomb_gauge_fixed(prec,func,fixed_conf,pars.gauge,pars.target_precision,&F_offset);
		get_out|=(not (iter<pars.nmax_iterations));
		get_out|=(not (iter%pars.unitarize_each==0));
		
		//switch off adaptative search if skipped too many times
		int nskipped_adapt_tol=5;
		if(use_adapt and nskipped_adapt>=nskipped_adapt_tol)
		  {
		    master_printf("Reached tolerance of skipping %d, switching off adaptative search\n",nskipped_adapt);
		    use_adapt=false;
		  }
		
		//switch off adaptative search if reached use_adapt_prec_tol
		double use_adapt_prec_tol=1e-14;
		if(use_adapt and prec<=use_adapt_prec_tol)
		  {
		    master_printf("Reached precision %lg, smaller than tolerance (%lg), switching off adaptative search\n",prec,use_adapt_prec_tol);
		    use_adapt=false;
		  }
	      }
	    while(not get_out);
	    
	    //now we put the fixer on su3, and make a real transformation
	    //on the basis of what we managed to fix
	    PAR(0,locVol,
		CAPTURE(TO_WRITE(fixer)),
		ivol,
		{
		  su3_unitarize_explicitly_inverting(fixer[ivol],fixer[ivol]);
		});
	    
	    gauge_transform_conf(fixed_conf,fixer,ext_conf);
	    
	    //check if really get out
	    really_get_out=false;
	    really_get_out|=check_Landau_or_Coulomb_gauge_fixed(prec,func,fixed_conf,pars.gauge,pars.target_precision,&F_offset);
	    really_get_out|=(not (iter<pars.nmax_iterations));
	  }
	while(not really_get_out);
	
	//crash if this did not work
	if(not really_get_out) crash("unable to fix to precision %16.16lg in %d iterations",pars.target_precision,pars.nmax_iterations);
	
	//free
	if(pars.use_generalized_cg) free_GCG_stuff();
	
	master_printf("Gauge fix time: %lg\n",time+take_time());
      }
  }
  
  void perform_random_gauge_transform(LxField<quad_su3>& conf_out,
				      const LxField<quad_su3>& conf_in)
  {
    LxField<su3> fixm("fixm",WITH_HALO);
    
    //extract random SU(3) matrix
    PAR(0,locVol,
	CAPTURE(b=maybeBackupLocRndGenForBenchmark(),
		TO_WRITE(fixm)),
	ivol,
	{
	  su3_put_to_rnd(fixm[ivol],loc_rnd_gen[ivol]);
	});
    
    //apply the transformation
    gauge_transform_conf(conf_out,fixm,conf_in);
  }
  
  void perform_random_gauge_transform(EoField<quad_su3>& conf_out,
				      const EoField<quad_su3>& conf_in)
  {
    
    //allocate fixing matrix
    EoField<su3> fixm("fixm",WITH_HALO);
    
    //extract random SU(3) matrix
    PAR(0,locVol,
	CAPTURE(b=maybeBackupLocRndGenForBenchmark(),
		TO_WRITE(fixm)),
	ivol,
	{
	  su3_put_to_rnd(fixm[loclx_parity[ivol]][loceo_of_loclx[ivol]],loc_rnd_gen[ivol]);
	});
    
    //apply the transformation
    gauge_transform_conf(conf_out,fixm,conf_in);
  }
}
