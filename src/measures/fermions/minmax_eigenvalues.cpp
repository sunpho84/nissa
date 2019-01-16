#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "dirac_operators/overlap/dirac_operator_overlap.hpp"
#include "eigenvalues/eigenvalues.hpp"
#include "geometry/geometry_mix.hpp"
#include "minmax_eigenvalues.hpp"
#include "new_types/rat_approx.hpp"
#include "operations/remez/remez_algorithm.hpp"

#include "dirac_operators/overlap/dirac_operator_overlap_kernel_portable.hpp"
#include "inverters/overlap/cgm_invert_overlap_kernel2.hpp"

#ifdef USE_THREADS
 #include "routines/thread.hpp"
#endif

///////////////////////////////////////////////
////      C. BONANNO AND M.CARDINALI       ////
///////////////////////////////////////////////

namespace nissa
{
  namespace minmax
  {
    THREADABLE_FUNCTION_4ARG(matrix_element_with_gamma, double*,out, complex*,buffer, spincolor*,x, int,igamma)
    {
      GET_THREAD_ID();
      
      NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
	{
	  spincolor t;
	  unsafe_dirac_prod_spincolor(t,base_gamma+igamma,x[ivol]);
	  spincolor_scalar_prod(buffer[ivol],x[ivol],t);
	}
      THREAD_BARRIER();
      
      complex_vector_glb_collapse(out,buffer,loc_vol);
    }
    THREADABLE_FUNCTION_END
  }
  
  //Computes the participation ratio
  double participation_ratio(spincolor *v)
  {
    GET_THREAD_ID();
    
    double *l=nissa_malloc("l",loc_vol,double);
    
    NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
      {
	complex t;
	spincolor_scalar_prod(t,v[ivol],v[ivol]);
	l[ivol]=t[RE];
      }
    THREAD_BARRIER();
    
    double s=double_vector_glb_norm2(l,loc_vol);
    double n2=double_vector_glb_norm2(v,loc_vol);
    
    return sqr(n2)/(glb_vol*s);
  }
  
  //measure minmax_eigenvalues
  void measure_minmax_eigenvalues(quad_su3 **conf_eo,theory_pars_t &theory_pars,minmax_eigenvalues_meas_pars_t &meas_pars,int iconf,int conf_created)
  {
    double eig_time=-take_time();
    
    FILE *fout=open_file(meas_pars.path,conf_created?"w":"a");
    
    //lx version
    quad_su3 *conf_lx=nissa_malloc("conf_lx",loc_vol+bord_vol,quad_su3);
    paste_eo_parts_into_lx_vector(conf_lx,conf_eo);
    
    rat_approx_t appr;
    int neigs=meas_pars.neigs;
    double residue=meas_pars.residue;
    double maxerr=sqrt(residue);
    
    //Parameters of the eigensolver
    const int mat_size=loc_vol*NCOL*NDIRAC;
    const int mat_size_to_allocate=(loc_vol+bord_vol)*NCOL*NDIRAC;
    const int niter_max=100000000;
    master_printf("mat_size=%d, mat_size_to_allocate=%d\n",mat_size,mat_size_to_allocate);
    
    //allocate
    complex *D_ov_eig_val=nissa_malloc("D_ov_eig_val",neigs,complex);
    complex **eigvec=nissa_malloc("eigvec",neigs,complex*);
    for(int ieig=0;ieig<neigs;ieig++)
      {
	eigvec[ieig]=nissa_malloc("eig",(loc_vol+bord_vol)*NDIRAC*NCOL,complex);
	vector_reset(eigvec[ieig]);
      }
    
    master_printf("neigs=%d, eig_precision=%.2e\n",neigs,maxerr);
    
    //consider only the first quark
    int iquark=0;
    double mass_overlap=theory_pars.quarks[iquark].mass_overlap;
    if(theory_pars.nflavs()!=1) crash("implemented only for 1 flavor");
    if(theory_pars.quarks[0].discretiz!=ferm_discretiz::OVERLAP) crash("Implemented only for overlap");
    
    rat_approx_for_overlap(conf_lx,&appr,mass_overlap,maxerr);
    
    //Verify the approximation
    {
      spincolor *in=(spincolor*)(eigvec[0]);
      generate_undiluted_source(in,RND_GAUSS,-1);
      spincolor *tmp=(spincolor*)(eigvec[1]);
      spincolor *out=(spincolor*)(eigvec[2]);
      summ_src_and_all_inv_overlap_kernel2_cgm(tmp,conf_lx,mass_overlap,&appr,niter_max,residue,in);
      apply_overlap_kernel(out,conf_lx,mass_overlap,tmp);
      summ_src_and_all_inv_overlap_kernel2_cgm(tmp,conf_lx,mass_overlap,&appr,niter_max,residue,out);
      apply_overlap_kernel(out,conf_lx,mass_overlap,tmp);
      
      double_vector_subtassign((double*)out,(double*)in,sizeof(spincolor)/sizeof(double)*loc_vol);
      
      double nout=double_vector_glb_norm2(out,loc_vol);
      double nin=double_vector_glb_norm2(in,loc_vol);
      
      master_printf("Norm of the source: %.16lg\n",sqrt(nin));
      master_printf("Norm of the difference: %.16lg\n",sqrt(nout));
      master_printf("Relative norm of the difference: %.16lg\n",sqrt(nout/nin));
    }
    
    appr.master_fprintf_expr(stdout);
    
    //Application of the Overlap Operator
    const auto imp_mat=[conf_lx,&theory_pars,&residue,iquark,&appr](complex *out_lx,complex *in_lx)
      {
	apply_overlap((spincolor*)out_lx,conf_lx,&appr,residue,theory_pars.quarks[iquark].mass_overlap,theory_pars.quarks[iquark].mass,(spincolor*)in_lx);
      };
    
    const auto filler=[](complex *out_lx){generate_undiluted_source((spincolor*)out_lx,RND_GAUSS,-1);};
    
    //Find eigenvalues and eigenvectors of the overlap
    eigenvalues_find(eigvec,D_ov_eig_val,neigs,meas_pars.min_max,mat_size,mat_size_to_allocate,imp_mat,maxerr,niter_max,filler);
    
    master_printf("\n\nEigenvalues of D Overlap:\n");
    for(int ieig=0;ieig<neigs;ieig++)
      for(FILE *f : {fout,stdout})
	master_fprintf(f,"%.16lg %.16lg\n",D_ov_eig_val[ieig][RE],D_ov_eig_val[ieig][IM]);
    
    complex *buffer=nissa_malloc("buffer",loc_vol,complex);
    spincolor *op=nissa_malloc("op",loc_vol,spincolor);
    spincolor *t=nissa_malloc("t",loc_vol+bord_vol,spincolor);
    
    //participation ratio
    master_printf("Participation ratio:\n");
    for(int ieig=0;ieig<neigs;ieig++)
      master_printf("%d: %.16lg\n",ieig,participation_ratio((spincolor*)(eigvec[ieig])));
    
    master_printf("Chirality of the eigenvectors:\n");
    for(int ieig=0;ieig<neigs;ieig++)
      {
	apply_overlap(op,conf_lx,&appr,residue,theory_pars.quarks[iquark].mass_overlap,theory_pars.quarks[iquark].mass,(spincolor*)(eigvec[ieig]));
	
	double n=double_vector_glb_norm2((spincolor*)(eigvec[ieig]),loc_vol);
	double e;
	double_vector_glb_scalar_prod(&e,(double*)(eigvec[ieig]),(double*)(op),sizeof(spincolor)/sizeof(double)*loc_vol);
	e/=n;
	
	//compute residue
	double_vector_summ_double_vector_prod_double((double*)op,(double*)op,(double*)(eigvec[ieig]),-e,sizeof(spincolor)/sizeof(double)*loc_vol);
	double r=sqrt(double_vector_glb_norm2(op,loc_vol)/n);
	
	//verify chirality
	complex c;
	minmax::matrix_element_with_gamma(c,buffer,(spincolor*)(eigvec[ieig]),5);
	
	//op=sign(H)*v
	spincolor *in=(spincolor*)(eigvec[ieig]);
	summ_src_and_all_inv_overlap_kernel2_cgm(t,conf_lx,mass_overlap,&appr,niter_max,residue,in);
	apply_overlap_kernel(op,conf_lx,mass_overlap,t);
	
	// i= (v, sign(H)*v) / |v|^2
	double i;
	double_vector_glb_scalar_prod(&i,(double*)(eigvec[ieig]),(double*)(op),sizeof(spincolor)/sizeof(double)*loc_vol);
	i/=n;
	
	// t=(1+g5 sign(H)) v
	GET_THREAD_ID();
	NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
	  {
	    unsafe_dirac_prod_spincolor(t[ivol],base_gamma+5,op[ivol]);
	    spincolor_summassign(t[ivol],op[ivol]);
	  }
	THREAD_BARRIER();
	
	// w= (v, (1+g5 sign(H))*v) / |v|^2
	double w;
	double_vector_glb_scalar_prod(&w,(double*)(eigvec[ieig]),(double*)(t),sizeof(spincolor)/sizeof(double)*loc_vol);
	w/=n;
	
	/// |(1+g5*sign(H))*v|^2
	double nop5=double_vector_glb_norm2(t,loc_vol);
	double h5=sqrt(nop5/n)-1;
	
	summ_src_and_all_inv_overlap_kernel2_cgm(t,conf_lx,mass_overlap,&appr,niter_max,residue,op);
	apply_overlap_kernel(op,conf_lx,mass_overlap,t);
	
	double_vector_subt((double*)t,(double*)op,(double*)in,sizeof(spincolor)/sizeof(double)*loc_vol);
	
	double nop=double_vector_glb_norm2(t,loc_vol);
	double h=sqrt(nop/n);
	
	master_printf(" %d (h: %lg, internal: %lg actual: %lg, hand: %lg, res: %lg, failure of sign(H): %lg, of g5*sign(H): %lg)\n   %lg %lg  %lg %lg %lg %lg\n",
		      ieig,i,D_ov_eig_val[ieig][RE],e,w,r,h,h5,c[RE],c[IM],((spincolor**)(eigvec))[ieig][0][0][0][RE],((spincolor**)eigvec)[ieig][0][1][0][RE],((spincolor**)eigvec)[ieig][0][2][0][RE],((spincolor**)eigvec)[ieig][0][3][0][RE]);
      }
    nissa_free(op);
    nissa_free(t);
    nissa_free(buffer);
    
    close_file(fout);
    
    eig_time+=take_time();
    master_printf("Eigenvalues time: %lg\n", eig_time);
    
    nissa_free(conf_lx);
    nissa_free(D_ov_eig_val);
    for(int ieig=0;ieig<neigs;ieig++) nissa_free(eigvec[ieig]);
    nissa_free(eigvec);
  }
  
  //print
  std::string minmax_eigenvalues_meas_pars_t::get_str(bool full)
  {
    std::ostringstream os;
    
    os<<"MeasMinMaxEigenval\n";
    os<<base_fermionic_meas_t::get_str(full);
    if(neigs!=def_neigs() or full) os<<" Neigs\t\t=\t"<<neigs<<"\n";
    if(wspace_size!=def_wspace_size() or full) os<<" WSpaceSize\t\t=\t"<<wspace_size<<"\n";
    if(min_max!=def_min_max() or full) os<<" MinMax\t\t=\t"<<min_max<<"\n";
    
    return os.str();
  }
}
