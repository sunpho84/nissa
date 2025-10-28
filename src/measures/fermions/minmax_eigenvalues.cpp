#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "dirac_operators/overlap/dirac_operator_overlap.hpp"
#include "eigenvalues/eigenvalues.hpp"
#include "eigenvalues/eigenvalues_overlap.hpp"
#include "geometry/geometry_mix.hpp"
#include "linalgs/reduce.hpp"
#include "minmax_eigenvalues.hpp"
#include "new_types/rat_approx.hpp"
#include "operations/remez/remez_algorithm.hpp"

#include "dirac_operators/overlap/dirac_operator_overlap_kernel_portable.hpp"
#include "inverters/overlap/cgm_invert_overlap_kernel2.hpp"

///////////////////////////////////////////////
////      C. BONANNO AND M.CARDINALI       ////
///////////////////////////////////////////////

namespace nissa
{
  namespace minmax
  {
    void matrix_element_with_gamma(double* out,complex* buffer,spincolor* x,int igamma)
    {
    CRASH("Reimplement");
      
    //   NISSA_PARALLEL_LOOP(ivol,0,locVol)
    // 	{
    // 	  spincolor t;
    // 	  unsafe_dirac_prod_spincolor(t,base_gamma[igamma],x[ivol]);
    // 	  spincolor_scalar_prod(buffer[ivol],x[ivol],t);
    // 	}
    //   NISSA_PARALLEL_LOOP_END;
    //   THREAD_BARRIER();
      
    //   glb_reduce((complex*)out,buffer,locVol);
    }
  }
  
  //Computes the participation ratio
  double participation_ratio(spincolor *v)
  {
    CRASH("Reimplement");
    
    // double *l=nissa_malloc("l",locVol,double);
    
    // NISSA_PARALLEL_LOOP(ivol,0,locVol)
    //   {
    // 	complex t;
    // 	spincolor_scalar_prod(t,v[ivol],v[ivol]);
    // 	l[ivol]=t[RE];
    //   }
    // NISSA_PARALLEL_LOOP_END;
    // THREAD_BARRIER();
    
    // double s=double_vector_glb_norm2(l,locVol);
    // double n2=double_vector_glb_norm2(v,locVol);
    
    // nissa_free(l);
    
    // return sqr(n2)/(glbVol*s);
  }
  
  //measure minmax_eigenvalues
  void measure_minmax_eigenvalues(eo_ptr<quad_su3> conf_eo,theory_pars_t &theory_pars,minmax_eigenvalues_meas_pars_t &meas_pars,int iconf,int conf_created)
  {
    CRASH("reimplement");
    // double eig_time=-take_time();
    
    // //Parameters of the eigensolver
    // FILE *fout=open_file(meas_pars.path,conf_created?"w":"a");
    // master_fprintf(fout," # iconf: %d\n",iconf);
    
    // //zero smooth time of the conf
    // quad_su3 *conf_lx=nissa_malloc("conf_lx",locVol+bord_vol,quad_su3);
    // // paste_eo_parts_into_lx_vector(conf_lx,conf_eo);
    // CRASH("reimplement");
    // //parameters of the measure
    // bool min_max=meas_pars.min_max;
    // int neigs=meas_pars.neigs;
    // double residue=meas_pars.residue;
    // int wspace_size=meas_pars.wspace_size;
    // double maxerr=sqrt(residue);
    
    // //allocate
    // complex *eigval=nissa_malloc("eigval",neigs,complex);
    // double *eig_res=nissa_malloc("eig_res",neigs,double);
    // spincolor **eigvec=nissa_malloc("eigvec",neigs,spincolor*);
    // spincolor *temp=nissa_malloc("temp",locVol,spincolor);
    // complex *buffer=nissa_malloc("buffer",locVol,complex);
    // for(int ieig=0;ieig<neigs;ieig++)
    //   eigvec[ieig]=nissa_malloc("eig",locVol+bord_vol,spincolor);
    
    // //loop on smooth
    // int nsmooth=0;
    // bool finished;
    // do
    //   {
    // 	VERBOSITY_LV1_MASTER_PRINTF("Measuring minmax_eigenvalues for nsmooth %d/%d\n",nsmooth,meas_pars.smooth_pars.nsmooth());
	
    // 	//plaquette for the current nsmooth
    // 	double plaq=global_plaquette_lx_conf(conf_lx);
    // 	master_fprintf(fout,"  # nsmooth: %d , plaq: %.16lg\n",nsmooth,plaq);
	
    // 	//loop on the quarks
    // 	for(int iquark=0;iquark<theory_pars.nflavs();iquark++)
    // 	  {
    // 	    master_fprintf(fout,"   # iquark: %d\n",iquark);
	    
    // 	    VERBOSITY_LV1_MASTER_PRINTF("Measuring minmax_eigenvalues for quark %d/%d\n",iquark+1,(int)theory_pars.nflavs());
	    
    // 	    double mass=theory_pars.quarks[iquark].mass;
    // 	    double mass_overlap=theory_pars.quarks[iquark].mass_overlap;
    // 	    if(theory_pars.quarks[0].discretiz!=ferm_discretiz::OVERLAP) CRASH("Implemented only for overlap");
	    
    // 	    //Generate the approximation
    // 	    rat_approx_t appr;
    // 	    generate_rat_approx_for_overlap(conf_lx,&appr,mass_overlap,maxerr);
    // 	    appr.master_fprintf_expr(stdout);
    // 	    verify_rat_approx_for_overlap(conf_lx,appr,mass_overlap,residue);
	    
    // 	    //Find the eigenvalues
    // 	    find_eigenvalues_overlap(eigvec,eigval,neigs,min_max,conf_lx,appr,residue,mass_overlap,mass,wspace_size);
	    
    // 	    //computes the participation ratio and chirality, recompute the eigenvalues and compute the residue
    // 	    for(int ieig=0;ieig<neigs;ieig++)
    // 	      {
    // 		master_fprintf(fout,"   # ieig: %d\n",ieig);
    // 		master_fprintf(fout,"     eigval: ( %.16lg , %.16lg )\n",eigval[ieig][RE],eigval[ieig][IM]);
		
    // 		//eigenvalue check
    // 		complex eigval_check;
    // 		apply_overlap(temp,conf_lx,&appr,residue,mass_overlap,mass,eigvec[ieig]);
    // 		double eigvec_norm2=double_vector_glb_norm2(eigvec[ieig],locVol);
    // 		complex_vector_glb_scalar_prod(eigval_check,(complex*)(eigvec[ieig]),(complex*)temp,sizeof(spincolor)/sizeof(complex)*locVol);
    // 		complex_prodassign_double(eigval_check,1.0/eigvec_norm2);
    // 		master_fprintf(fout,"     eigval_check: ( %.16lg , %.16lg)\n",eigval_check[RE],eigval_check[IM]);
		
    // 		//residue
    // 		complex_vector_subtassign_complex_vector_prod_complex((complex*)temp,(complex*)(eigvec[ieig]),eigval_check,sizeof(spincolor)/sizeof(complex)*locVol);
    // 		eig_res[ieig]=sqrt(double_vector_glb_norm2(temp,locVol)/eigvec_norm2);
    // 		master_fprintf(fout,"     residue: %.16lg\n",residue);
		
    // 		//participation ratio
    // 		double pr=participation_ratio(eigvec[ieig]);
    // 		master_fprintf(fout,"     partic_rat: %.16lg\n",pr);
		
    // 		//chirality
    // 		complex chir;
    // 		minmax::matrix_element_with_gamma(chir,buffer,eigvec[ieig],5);
    // 		master_fprintf(fout,"     chirality: %.16lg %.16lg\n",chir[RE],chir[IM]);
		
    // 		MASTER_PRINTF("\n");
    // 	      }
    // 	  }
	
    // 	//proceeds with smoothing
    // 	CRASH("reimplement");	  finished=1;
    // 	// finished=smooth_lx_conf_until_next_meas(conf_lx,meas_pars.smooth_pars,nsmooth);
    //   }
    // while(not finished);
    
    // //close the file
    // close_file(fout);
    
    // //print elapsed time
    // eig_time+=take_time();
    // MASTER_PRINTF("Eigenvalues computation time: %lg\n", eig_time);
    
    // //free
    // nissa_free(conf_lx);
    // nissa_free(eigval);
    // nissa_free(eig_res);
    // nissa_free(temp);
    // nissa_free(buffer);
    // for(int ieig=0;ieig<neigs;ieig++)
    //   nissa_free(eigvec[ieig]);
    // nissa_free(eigvec);
  }
}
