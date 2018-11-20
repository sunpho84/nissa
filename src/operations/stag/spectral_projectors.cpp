#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <algorithm>

#include "base/thread_macros.hpp"
#include "base/vectors.hpp"
#include "communicate/borders.hpp"
#include "communicate/edges.hpp"
#include "geometry/geometry_mix.hpp"
#include "io/input.hpp"
#include "linalgs/linalgs.hpp"
#include "new_types/complex.hpp"
#include "new_types/float_128.hpp"
#include "new_types/spin.hpp"
#include "new_types/su3.hpp"
#include "operations/fft.hpp"
#include "operations/gaugeconf.hpp"
#include "operations/remap_vector.hpp"
#include "operations/su3_paths/gauge_sweeper.hpp"
#include "operations/smearing/stout.hpp"
#include "operations/smearing/Wflow.hpp"
#include "operations/su3_paths/plaquette.hpp"
#include "routines/ios.hpp"
#include "routines/mpi_routines.hpp"
#include "dirac_operators/tmQ/dirac_operator_tmQ.hpp"
#include "dirac_operators/tmQ2/dirac_operator_tmQ2.hpp"
#include "dirac_operators/stD/dirac_operator_stD.hpp"
#include "eigenvalues/eigenvalues_all.hpp"
#include "eigenvalues/eigenvalues.hpp"
#include "eigenvalues/eigenvalues_autarchic.hpp"
#include "inverters/twisted_mass/cg_invert_tmQ2.hpp"

//FIXME clean includes
#include "base/random.hpp"
#include "base/vectors.hpp"
#include "operations/stag/stag.hpp"
#include "communicate/borders.hpp"
#include "geometry/geometry_eo.hpp"
#include "hmc/backfield.hpp"
#include "inverters/staggered/cg_invert_stD.hpp"
#include "linalgs/linalgs.hpp"
#include "new_types/su3.hpp"
#include "routines/mpi_routines.hpp"

#ifdef USE_THREADS
 #include "routines/thread.hpp"
#endif

#include "spectral_projectors.hpp"

namespace nissa{
  // This measure will compute the first 'n' eigenvalues (parameter) 
  // and eigenvectors of DD^+ in staggered formulation, in order to build an estimate
  // of the topological charge as Q = \sum_i \bar(u)_i gamma5 u_i,
  // where {u_i} are the first n eigenvectors.

  THREADABLE_FUNCTION_4ARG(fill_eigenpart, complex**,eigvec, quad_su3**,conf,int,neigs, double,eig_precision)
  {
//    GET_THREAD_ID();
    master_fprintf(stdout,"neigs = %d, eig_precision = %.2e\n", neigs, eig_precision);

//    complex lambda[neig];

//    //wrap the application of DD^+ into an object that can be passed to the eigenfinder
    color *out_tmp_eo[2]={nissa_malloc("out_tmp_EVN",loc_volh+bord_volh,color),nissa_malloc("out_tmp_ODD",loc_volh+bord_volh,color)};
    color *in_tmp_eo[2]={nissa_malloc("in_tmp_EVN",loc_volh+bord_volh,color),nissa_malloc("in_tmp_ODD",loc_volh+bord_volh,color)};
    color *fill_tmp_eo[2]={nissa_malloc("fill_tmp_EVN",loc_volh+bord_volh,color),nissa_malloc("fill_tmp_ODD",loc_volh+bord_volh,color)};
		
		if(bord_volh!=0)
			master_fprintf(stderr,"problem! bord_volh!=0 unimplemented\n");

		const auto imp_mat=[conf,in_tmp_eo,out_tmp_eo](complex *out,complex *in){
        // apply_tmQ2((spincolor*)out,conf,kappa,temp_imp_mat,mu,(spincolor*)in);
        
        // apply_tmQ((spincolor*)out,conf,kappa,mu,(spincolor*)in);
      GET_THREAD_ID();
      NISSA_PARALLEL_LOOP(ivol,0,loc_volh)
      {
			  in_tmp_eo[EVN][ivol][0][RE]=in[ivol*3][RE];
			  in_tmp_eo[EVN][ivol][0][IM]=in[ivol*3][IM];
			  in_tmp_eo[EVN][ivol][1][RE]=in[ivol*3+1][RE];
			  in_tmp_eo[EVN][ivol][1][IM]=in[ivol*3+1][IM];
			  in_tmp_eo[EVN][ivol][2][RE]=in[ivol*3+2][RE];
			  in_tmp_eo[EVN][ivol][2][IM]=in[ivol*3+2][IM];
      
			  in_tmp_eo[ODD][ivol][0][RE]=in[(ivol+loc_volh)*3][RE];
			  in_tmp_eo[ODD][ivol][0][IM]=in[(ivol+loc_volh)*3][IM];
			  in_tmp_eo[ODD][ivol][1][RE]=in[(ivol+loc_volh)*3+1][RE];
			  in_tmp_eo[ODD][ivol][1][IM]=in[(ivol+loc_volh)*3+1][IM];
			  in_tmp_eo[ODD][ivol][2][RE]=in[(ivol+loc_volh)*3+2][RE];
			  in_tmp_eo[ODD][ivol][2][IM]=in[(ivol+loc_volh)*3+2][IM];
      }
						
      evn_apply_stD(out_tmp_eo[EVN],conf,0.0,(color**)in_tmp_eo);		
      odd_apply_stD(out_tmp_eo[ODD],conf,0.0,(color**)in_tmp_eo);
//      std::swap(out_tmp_eo,in_tmp_eo);    
      evn_apply_stD_dag(in_tmp_eo[EVN],conf,0.0,(color**)out_tmp_eo);		
      odd_apply_stD_dag(in_tmp_eo[ODD],conf,0.0,(color**)out_tmp_eo);
//        apply_stD(out_tmp_eo,conf,mu,in_tmp_eo);
//
//        NISSA_PARALLEL_LOOP(ivol,0,loc_volh)
//					out[ivol]=out_tmp_eo[EVN][ivol];
//        NISSA_PARALLEL_LOOP(ivol,0,loc_volh)
//					out[ivol+loc_volh]=out_tmp_eo[ODD][ivol];
      NISSA_PARALLEL_LOOP(ivol,0,loc_volh)
      {
			  out[ivol*3  ][RE] =in_tmp_eo[EVN][ivol][0][RE];
			  out[ivol*3  ][IM] =in_tmp_eo[EVN][ivol][0][IM];
			  out[ivol*3+1][RE] =in_tmp_eo[EVN][ivol][1][RE];
			  out[ivol*3+1][IM] =in_tmp_eo[EVN][ivol][1][IM];
			  out[ivol*3+2][RE] =in_tmp_eo[EVN][ivol][2][RE];
			  out[ivol*3+2][IM] =in_tmp_eo[EVN][ivol][2][IM];
			  
        out[(ivol+loc_volh)*3  ][RE] =in_tmp_eo[ODD][ivol][0][RE];
			  out[(ivol+loc_volh)*3  ][IM] =in_tmp_eo[ODD][ivol][0][IM];
			  out[(ivol+loc_volh)*3+1][RE] =in_tmp_eo[ODD][ivol][1][RE];
			  out[(ivol+loc_volh)*3+1][IM] =in_tmp_eo[ODD][ivol][1][IM];
			  out[(ivol+loc_volh)*3+2][RE] =in_tmp_eo[ODD][ivol][2][RE];
			  out[(ivol+loc_volh)*3+2][IM] =in_tmp_eo[ODD][ivol][2][IM];
      }
				
//
//        GET_THREAD_ID();
//        NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
//          unsafe_dirac_prod_spincolor(((color*)out)[ivol],base_gamma+5,temp_imp_mat[ivol]);
        set_borders_invalid(out);
      };
  


////    const auto imp_mat=[conf,mu=am,temp_imp_mat](complex *out,complex *in)
////      {
////        apply_stD(temp_imp_mat,conf,mu,(color*)in);
////        GET_THREAD_ID();
////        NISSA_PARALLEL_LOOP(ivol,0,loc_vol)
////          unsafe_dirac_prod_spincolor(((spincolor*)out)[ivol],base_gamma+5,temp_imp_mat[ivol]);
////        set_borders_invalid(out);
////      };
//    
    //parameters of the eigensolver
    const bool min_max=0;
    const int mat_size=loc_vol*sizeof(color)/sizeof(complex);
    const int mat_size_to_allocate=(loc_vol+bord_vol)*sizeof(color)/sizeof(complex);
    const int niter_max=50; //XXX changed from 100000
    master_printf("mat_size=%d, mat_size_to_allocate=%d\n",mat_size,mat_size_to_allocate); 
    //wrap the generation of the test vector into an object that can be passed to the eigenfinder
    const auto filler=[fill_tmp_eo](complex *a)
      {
        generate_fully_undiluted_eo_source((color**)fill_tmp_eo,RND_GAUSS,-1,0); 
        GET_THREAD_ID();
        NISSA_PARALLEL_LOOP(ivol,0,loc_volh)
        {
          a[ivol*3  ][RE] =fill_tmp_eo[EVN][ivol][0][RE];
          a[ivol*3  ][IM] =fill_tmp_eo[EVN][ivol][0][IM];
          a[ivol*3+1][RE] =fill_tmp_eo[EVN][ivol][1][RE];
          a[ivol*3+1][IM] =fill_tmp_eo[EVN][ivol][1][IM];
          a[ivol*3+2][RE] =fill_tmp_eo[EVN][ivol][2][RE];
          a[ivol*3+2][IM] =fill_tmp_eo[EVN][ivol][2][IM];
          
          a[(ivol+loc_volh)*3  ][RE] =fill_tmp_eo[ODD][ivol][0][RE];
          a[(ivol+loc_volh)*3  ][IM] =fill_tmp_eo[ODD][ivol][0][IM];
          a[(ivol+loc_volh)*3+1][RE] =fill_tmp_eo[ODD][ivol][1][RE];
          a[(ivol+loc_volh)*3+1][IM] =fill_tmp_eo[ODD][ivol][1][IM];
          a[(ivol+loc_volh)*3+2][RE] =fill_tmp_eo[ODD][ivol][2][RE];
          a[(ivol+loc_volh)*3+2][IM] =fill_tmp_eo[ODD][ivol][2][IM];
        }
      };
    
    //launch the eigenfinder
    double eig_time=-take_time();
    complex DD_eig_val[neigs];
    eigenvalues_of_hermatr_find(eigvec,DD_eig_val,neigs,min_max,mat_size,mat_size_to_allocate,imp_mat,eig_precision,niter_max,filler);

    double charge_cut[neigs];

    master_printf("\n\nEigenvalues of Id:\n");
    for(int ieig=0; ieig<neigs; ++ieig){
      master_printf("%d (%.16lg,%.16lg)\n",ieig,DD_eig_val[ieig][RE],DD_eig_val[ieig][IM]);
      
      // compute partial sum of tr(g5)      
      // save 'eigvec[ieigs]' in staggered format (i.e., 'fill_tmp_eo'),
      // then multiply it with gamma5 and save the result in 
      // 'in_tmp_eo'. The term contributing to the partial sum of tr(g5)
      // will be the hermitian product between 'fill_tmp_eo' and 'in_tmp_eo'  

      // 'complex*' to 'color**' conversion (TODO: industrialize)
      GET_THREAD_ID();
      NISSA_PARALLEL_LOOP(ivol,0,loc_volh)
      {
        fill_tmp_eo[EVN][ivol][0][RE]=eigvec[ieig][ivol*3][RE];
        fill_tmp_eo[EVN][ivol][0][IM]=eigvec[ieig][ivol*3][IM];
        fill_tmp_eo[EVN][ivol][1][RE]=eigvec[ieig][ivol*3+1][RE];
        fill_tmp_eo[EVN][ivol][1][IM]=eigvec[ieig][ivol*3+1][IM];
        fill_tmp_eo[EVN][ivol][2][RE]=eigvec[ieig][ivol*3+2][RE];
        fill_tmp_eo[EVN][ivol][2][IM]=eigvec[ieig][ivol*3+2][IM];

        fill_tmp_eo[ODD][ivol][0][RE]=eigvec[ieig][(ivol+loc_volh)*3][RE];
        fill_tmp_eo[ODD][ivol][0][IM]=eigvec[ieig][(ivol+loc_volh)*3][IM];
        fill_tmp_eo[ODD][ivol][1][RE]=eigvec[ieig][(ivol+loc_volh)*3+1][RE];
        fill_tmp_eo[ODD][ivol][1][IM]=eigvec[ieig][(ivol+loc_volh)*3+1][IM];
        fill_tmp_eo[ODD][ivol][2][RE]=eigvec[ieig][(ivol+loc_volh)*3+2][RE];
        fill_tmp_eo[ODD][ivol][2][IM]=eigvec[ieig][(ivol+loc_volh)*3+2][IM];
      }

      //multiply by gamma5
//      apply_stag_op(out_tmp_eo,conf,u1b,GAMMA_5,IDENTITY,fill_tmp_eo);

      //take hermitian product (TODO) 
    }
    master_printf("\n\n\n");

    eig_time+=take_time();
    master_printf("Eigenvalues time: %lg\n",eig_time);
    
    nissa_free(out_tmp_eo[EVN]);
    nissa_free(out_tmp_eo[ODD]);
    nissa_free(in_tmp_eo[EVN]);
    nissa_free(in_tmp_eo[ODD]);
    nissa_free(fill_tmp_eo[EVN]);
    nissa_free(fill_tmp_eo[ODD]);
  }
  THREADABLE_FUNCTION_END

  
  //measure the topological charge
  void measure_spectral_proj(quad_su3 **conf,theory_pars_t &theory_pars,spectr_proj_meas_pars_t &meas_pars, int iconf,bool conf_created)
  {
    FILE *file=open_file(meas_pars.path,conf_created?"w":"a");
//    FILE *file_proj=open_file(meas_pars.path+"%s_proj_x",conf_created?"w":"a");
//    start_glb_rnd_gen(1020419);
//    start_loc_rnd_gen(1292024);
    int neigs = meas_pars.neigs;
//    int ncopies=meas_pars.ncopies;
    
    
//    double *charge=nissa_malloc("charge",loc_vol,double);
    

    complex *eigvec[neigs];
    master_printf("loc_vol*3=%d, ",loc_vol*3);
    for(int ieig=0;ieig<neigs;ieig++){
		  eigvec[ieig]=nissa_malloc("eigvec",loc_vol*3,complex);
    }
//    for(int icopy=0;icopy<ncopies;icopy++)
//      {
//        master_fprintf(file,"%d",iconf);
         
        //measure magnetization for each quark
//        for(int iflav=0;iflav<theory_pars.nflavs();iflav++)
//          {
//          if(theory_pars.quarks[iflav].discretiz!=ferm_discretiz::ROOT_STAG) crash("not defined for non-staggered quarks");
	    
//            complex magn={0,0};
//            complex magn_proj_x[glb_size[1]]; //this makes pair and pact with "1" and "2" upstairs
//            for(int i=0;i<glb_size[1];i++) magn_proj_x[i][RE]=magn_proj_x[i][IM]=0;
            
            //loop over hits
            int nhits=meas_pars.nhits;
            for(int hit=0;hit<nhits;hit++)
              {
                verbosity_lv2_master_printf("Evaluating spectral topological charge nhits %d/%d\n",hit+1,nhits);
            
                //compute and summ
                //complex temp,temp_magn_proj_x[glb_size[1]];
                
                

                fill_eigenpart(eigvec,conf,meas_pars.neigs,meas_pars.eig_precision);
                

//                magnetization(&temp,temp_magn_proj_x,meas_pars.rnd_type,conf,theory_pars.em_field_pars.flag,theory_pars.backfield[iflav],&theory_pars.quarks[iflav],meas_pars.residue); //flag holds quantization
//                
//                //normalize
//                complex_summ_the_prod_double(magn,temp,1.0/nhits);
//                for(int x=0;x<glb_size[1];x++) complex_summ_the_prod_double(magn_proj_x[x],temp_magn_proj_x[x],1.0/nhits);
              }
            
            //output
//            master_fprintf(file,"\t%+016.16lg \t%+016.16lg",magn[RE],magn[IM]);
//            master_fprintf(file,"ok\n");
//            for(int x=0;x<glb_size[1];x++)
//              master_fprintf(file_proj,"%d\t%d\t%d\t%d\t%+016.16lg \t%+016.16lg\n",iconf,icopy,iflav,x,magn_proj_x[x][RE],magn_proj_x[x][IM]);
//          }
//        
//        master_fprintf(file,"\n");
//      }
    
//    close_file(file_proj);
  
  //print

//      fill_eigenpart(eigvec_conv,eigvec,neig,unsmoothed_conf,pars.kappa,pars.am,pars.r,pars.eig_precision);
//      finished=true;
//      {}while(!finished);
    
    //discard smoothed conf
//    nissa_free(charge);
//   	
//    for(int ieig=0;ieig<neigs;ieig++){
//			nissa_free(eigvec[ieig][EVN]);
//			nissa_free(eigvec[ieig][ODD]);
//      nissa_free(eigvec[ieig]);
////			nissa_free(eigvec_conv[EVN][ieig]);
////			nissa_free(eigvec_conv[ODD][ieig]);
//	  }
    
    for(int ieig=0; ieig<neigs; ieig++){
      nissa_free(eigvec[ieig]);
    }

    close_file(file);
  }
  
  //print pars
  std::string spectr_proj_meas_pars_t::get_str(bool full)
  {
    std::ostringstream os;
    
    os<<"MeasSpectralProj\n";
    os<<base_fermionic_meas_t::get_str(full);
    if(neigs!=def_neigs() or full) os<<" Neigs\t\t=\t"<<neigs<<"\n";
    if(eig_precision!=def_eig_precision() or full) os<<" EigPrecision\t\t=\t"<<eig_precision<<"\n";
    
    return os.str();
  }
}
