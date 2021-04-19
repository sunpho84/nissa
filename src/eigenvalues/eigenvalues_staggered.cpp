#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <random/randomGenerate.hpp>
#include <dirac_operators/stD/dirac_operator_stD.hpp>
#include <eigenvalues/eigenvalues.hpp>
#include <geometry/geometry_eo.hpp>
#include <geometry/geometry_mix.hpp>
#include <hmc/backfield.hpp>

namespace nissa
{
  //computes the spectrum of the staggered operator
  void find_eigenvalues_staggered_D2ee(color **eigvec,complex *eigval,int neigs,bool min_max,eo_ptr<quad_su3> conf,eo_ptr<quad_u1> u1b,double mass2,double residue,int wspace_size)
  {
    add_backfield_with_stagphases_to_conf(conf,u1b);
    
    color *temp=nissa_malloc("temp",locVolhWithBord.nastyConvert(),color);
    
    //Application of the staggered Operator
    const auto imp_mat=[conf,&temp,mass2](complex *out_e,complex *in_e)
      {
	apply_stD2ee_m2((color*)out_e,conf,temp,mass2,(color*)in_e);
      };
    
    const auto filler=[](complex *out_e)
      {
	generate_fully_undiluted_eo_source((color*)out_e,RND_GAUSS,-1,EVN);
      };
    
    //parameters of the finder
    const int mat_size=locVolh()*NCOL;
    const int mat_size_to_allocate=locVolhWithBord()*NCOL;
    const int niter_max=100000000;
    master_printf("mat_size=%d, mat_size_to_allocate=%d\n",mat_size,mat_size_to_allocate);
    
    //precision of the eigenvalues
    double maxerr=sqrt(residue);
    
    verbosity_lv1_master_printf("Starting to search for %d %s eigenvalues of the Staggered operator DD^+ even-projected, with a precision of %lg, and Krylov space size of %d\n",neigs,(min_max?"max":"min"),maxerr,wspace_size);
    
    //find eigenvalues and eigenvectors of the staggered
    eigenvalues_find((complex**)eigvec,eigval,neigs,min_max,mat_size,mat_size_to_allocate,imp_mat,maxerr,niter_max,filler,wspace_size);
    
    rem_backfield_with_stagphases_from_conf(conf,u1b);
    
    nissa_free(temp);
  }
  
  //computes the spectrum of the staggered iD operator
  void find_eigenvalues_staggered_iD(color **eigvec,complex *eigval,int neigs,bool min_max,eo_ptr<quad_su3> conf,eo_ptr<quad_u1> u1b,double residue,int wspace_size)
  {
    eo_ptr<color> temp_in_eo={nissa_malloc("temp_in_eo_EVN",locVolhWithBord.nastyConvert(),color),nissa_malloc("temp_in_eo_ODD",locVolhWithBord.nastyConvert(),color)};
    eo_ptr<color> temp_out_eo={nissa_malloc("temp_out_eo_EVN",locVolhWithBord.nastyConvert(),color),nissa_malloc("temp_out_eo_ODD",locVolhWithBord.nastyConvert(),color)};
    
    const auto imp_mat=[conf,u1b,&temp_in_eo,&temp_out_eo](complex *out,complex *in)
      {
        split_lx_vector_into_eo_parts(temp_in_eo,(color*)in);
        
        // temp_out_eo = D * in
        add_backfield_with_stagphases_to_conf(conf,u1b);
        apply_stD(temp_out_eo,conf,0.0,temp_in_eo);
        rem_backfield_with_stagphases_from_conf(conf,u1b);
	
	for(int eo=0;eo<2;eo++)
	  {
	    
	    // temp_out_eo = i * D * in
	    NISSA_PARALLEL_LOOP(ieo,0,locVolh)
	      for(int ic=0;ic<NCOL;ic++)
		assign_complex_prod_i(temp_out_eo[eo][ieo.nastyConvert()][ic]);
	    NISSA_PARALLEL_LOOP_END;
	  }
	
        paste_eo_parts_into_lx_vector((color*)out,temp_out_eo);
      };
    
    const auto filler=
      [](complex *out)
      {
	generate_fully_undiluted_eo_source((color*)out,RND_GAUSS,-1,EVN);
      };
    
    //parameters of the finder
    const int mat_size=2*locVolh()*NCOL;
    const int mat_size_to_allocate=2*(locVolhWithBord.nastyConvert())*NCOL;
    const int niter_max=100000000;
    master_printf("mat_size=%d, mat_size_to_allocate=%d\n",mat_size,mat_size_to_allocate);
    
    //precision of the eigenvalues
    double maxerr=sqrt(residue);
    
    verbosity_lv1_master_printf("Starting to search for %d %s eigenvalues of the Staggered iD operator with a precision of %lg, and Krylov space size of %d\n",
				neigs,(min_max?"max":"min"),maxerr,wspace_size);
    
    //find eigenvalues and eigenvectors of the staggered
    eigenvalues_find((complex**)eigvec,eigval,neigs,min_max,mat_size,mat_size_to_allocate,imp_mat,maxerr,niter_max,filler,wspace_size);
    
    nissa_free(temp_in_eo[EVN]);
    nissa_free(temp_in_eo[ODD]);
    nissa_free(temp_out_eo[EVN]);
    nissa_free(temp_out_eo[ODD]);
  }
  
  //computes the spectrum of the staggered Adams operator (iD_st - Gamma5 m_Adams)
  void find_eigenvalues_staggered_Adams(color **eigvec,complex *eigval,int neigs,bool min_max,eo_ptr<quad_su3> conf,eo_ptr<quad_u1> u1b,double mass,double m_Adams,double residue,int wspace_size)
  {
    eo_ptr<color> temp={nissa_malloc("temp_EVN",locVolhWithBord.nastyConvert(),color),nissa_malloc("temp_ODD",locVolhWithBord.nastyConvert(),color)};
    eo_ptr<color> temp_in_eo = {nissa_malloc("temp_in_eo_EVN",locVolhWithBord.nastyConvert(),color),nissa_malloc("temp_in_eo_ODD",locVolhWithBord.nastyConvert(),color)};
    eo_ptr<color> temp_out_eo = {nissa_malloc("temp_out_eo_EVN",locVolhWithBord.nastyConvert(),color),nissa_malloc("temp_out_eo_ODD",locVolhWithBord.nastyConvert(),color)};
    
    //Application of the staggered Operator
    const auto imp_mat=[conf,u1b,&temp,&temp_in_eo,&temp_out_eo,mass,m_Adams](complex *out,complex *in)
      {
        split_lx_vector_into_eo_parts(temp_in_eo,(color*)in);
        
        apply_Adams(temp_out_eo,conf,u1b,mass,m_Adams,temp,temp_in_eo);
	
        paste_eo_parts_into_lx_vector((color*)out,temp_out_eo);
      };
    
    const auto filler=
      [](complex *out)
      {
	generate_fully_undiluted_eo_source((color*)out,RND_GAUSS,-1,EVN);
      };
    
    //parameters of the finder
    const int mat_size=2*locVolh()*NCOL;
    const int mat_size_to_allocate=2*(locVolhWithBord.nastyConvert())*NCOL;
    const int niter_max=100000000;
    master_printf("mat_size=%d, mat_size_to_allocate=%d\n",mat_size,mat_size_to_allocate);
    
    //precision of the eigenvalues
    double maxerr=sqrt(residue);
    
    verbosity_lv1_master_printf("Starting to search for %d %s eigenvalues of the Staggered Adams operator with a precision of %lg, and Krylov space size of %d\n",neigs,(min_max?"max":"min"),maxerr,wspace_size);
    
    //find eigenvalues and eigenvectors of the staggered
    eigenvalues_find((complex**)eigvec,eigval,neigs,min_max,mat_size,mat_size_to_allocate,imp_mat,maxerr,niter_max,filler,wspace_size);
    
    nissa_free(temp[EVN]);
    nissa_free(temp[ODD]);
    nissa_free(temp_in_eo[EVN]);
    nissa_free(temp_in_eo[ODD]);
    nissa_free(temp_out_eo[EVN]);
    nissa_free(temp_out_eo[ODD]);
  }
  
  //computes the spectrum of the staggered Adams operator (Eps D_st - Gamma5 m_Adams), where Eps = Gamma5 x Gamma5.
  void find_eigenvalues_staggered_AdamsII(color **eigvec,complex *eigval,int neigs,bool min_max,eo_ptr<quad_su3> conf,eo_ptr<quad_u1> u1b,double mass,double m_Adams,double residue,int wspace_size)
  {
    eo_ptr<color> temp={nissa_malloc("temp_EVN",locVolhWithBord.nastyConvert(),color),nissa_malloc("temp_ODD",locVolhWithBord.nastyConvert(),color)};
    eo_ptr<color> temp_in_eo={nissa_malloc("temp_in_eo_EVN",locVolhWithBord.nastyConvert(),color),nissa_malloc("temp_in_eo_ODD",locVolhWithBord.nastyConvert(),color)};
    eo_ptr<color> temp_out_eo={nissa_malloc("temp_out_eo_EVN",locVolhWithBord.nastyConvert(),color),nissa_malloc("temp_out_eo_ODD",locVolhWithBord.nastyConvert(),color)};
    
    //Application of the staggered Operator
    const auto imp_mat=
      [conf,u1b,&temp,&temp_in_eo,&temp_out_eo,mass,m_Adams](complex *out,complex *in)
      {
        split_lx_vector_into_eo_parts(temp_in_eo,(color*)in);
        
        apply_AdamsII(temp_out_eo,conf,u1b,mass,m_Adams,temp,temp_in_eo);
	
        paste_eo_parts_into_lx_vector((color*)out,temp_out_eo);
      };
    
    const auto filler=
      [](complex *out)
      {
	generate_fully_undiluted_eo_source((color*)out,RND_GAUSS,-1,EVN);
      };
    
    //parameters of the finder
    const int mat_size=2*locVolh()*NCOL;
    const int mat_size_to_allocate=2*(locVolhWithBord.nastyConvert())*NCOL;
    const int niter_max=100000000;
    master_printf("mat_size=%d, mat_size_to_allocate=%d\n",mat_size,mat_size_to_allocate);
    
    //precision of the eigenvalues
    double maxerr=sqrt(residue);
    
    verbosity_lv1_master_printf("Starting to search for %d %s eigenvalues of the Staggered Adams operator with a precision of %lg, and Krylov space size of %d\n",neigs,(min_max?"max":"min"),maxerr,wspace_size);
    
    //find eigenvalues and eigenvectors of the staggered
    eigenvalues_find((complex**)eigvec,eigval,neigs,min_max,mat_size,mat_size_to_allocate,imp_mat,maxerr,niter_max,filler,wspace_size);
    
    nissa_free(temp[EVN]);
    nissa_free(temp[ODD]);
    nissa_free(temp_in_eo[EVN]);
    nissa_free(temp_in_eo[ODD]);
    nissa_free(temp_out_eo[EVN]);
    nissa_free(temp_out_eo[ODD]);
  }
}
