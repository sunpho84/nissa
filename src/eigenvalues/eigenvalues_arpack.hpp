#ifndef _EIGENVALUES_ARPACK_HPP
#define _EIGENVALUES_ARPACK_HPP

#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <mpi.h>

#include <arpack/parpack.hpp>

#include "base/vectors.hpp"
#include "new_types/complex.hpp"
#include "routines/ios.hpp"

#ifndef EXTERN_ARPACK
 #define EXTERN_ARPACK extern
 #define INIT_ARPACK_TO(...)
#else
 #define INIT_ARPACK_TO(...) __VA_ARGS__
#endif

namespace nissa
{
  namespace arpack_data
  {
    //normal eigenvalue problem
    EXTERN_ARPACK char bmat[] INIT_ARPACK_TO(="I");
    //largest or smallest
    EXTERN_ARPACK char which[2][3] INIT_ARPACK_TO(={"SM","LM"});
    //parameters
    EXTERN_ARPACK int iparam[11];
    //pointer to work
    EXTERN_ARPACK int ipntr[14];
    
    //type for cast
    typedef _Complex double c_t;
    
    void check_exit_status_arpack(const int info);
    void init_params(int *iparam,int niter_max);
    void eigenvalues_of_hermatr_find_arpack_compute(MPI_Fint &comm,complex **eig_vec,double *eig_val,int neig,bool min_max,
						    const int mat_size,const double target_precision,const int wspace_size,
						    complex *residue,complex *workd,complex *workl,const int lworkl,complex *rwork,int &info,
						    complex *v,const int ldv);
  }
  
  //use arpack
  template <class Fmat,class Filler>
  void eigenvalues_of_hermatr_find_arpack(complex **eig_vec,double *eig_val,int neig,bool min_max,
					  const int mat_size,const int mat_size_to_allocate,const Fmat &imp_mat,
					  const double target_precision,const int niter_max,
					  const Filler &filler)
  {
    using namespace arpack_data;
    
    //convert the communicator
    MPI_Fint comm=MPI_Comm_c2f(MPI_COMM_WORLD);
    
    //vector of residue: we fill it
    complex *residue=nissa_malloc("residue",mat_size_to_allocate,complex);
    // filler(residue);
    // int info=1; //do not use arpack random generator
    int info=0; //use arpack random generator
    
    //base for eigenvectors
    const int wspace_size=2*neig;
    //workspace
    complex *v=nissa_malloc("v",wspace_size*mat_size_to_allocate,complex);
    //size of each line
    const int ldv=mat_size_to_allocate;
    //internal workspace
    complex *workd=nissa_malloc("workd",3*mat_size,complex);
    const int lworkl=(3*wspace_size+5)*wspace_size;
    complex *workl=nissa_malloc("workl",lworkl,complex);
    complex *rwork=nissa_malloc("rwork",wspace_size,complex);
    
    //initialize the parameters
    init_params(iparam,niter_max);
    
    //temporary vectors
    complex *temp_x=nissa_malloc("temp_x",mat_size_to_allocate,complex);
    complex *temp_y=nissa_malloc("temp_y",mat_size_to_allocate,complex);
    
    //main loop to find Ritz basis
    bool goon=true;
    int iter=0;
    do
      {
	//invoke the step
	int ido;
	arpack::internal::pznaupd_c(comm,&ido,bmat,mat_size,which[min_max],neig,target_precision,(c_t*)residue,wspace_size,(c_t*)v,ldv,iparam,ipntr,(c_t*)workd,(c_t*)workl,lworkl,(c_t*)rwork,&info);
	verbosity_lv1_master_printf("iteration %d, ido: %d, info: %d\n",iter,ido,info);
	
	switch(ido)
	  {
	    complex *x,*y;
	  case 1:
	    goon=true;
	    
	    //reverse communication
	    x=workd+ipntr[0]-1;
	    y=workd+ipntr[1]-1;
	    
	    memcpy(temp_x,x,sizeof(complex)*mat_size);
	    set_borders_invalid(temp_x);
	    imp_mat(temp_y,temp_x);
	    memcpy(y,temp_y,sizeof(complex)*mat_size);
	    break;
	  case 99:
	    goon=false;
	    master_printf("Finished!\n");
	  default:
	    crash("ido %d not contemplated",ido);
	  }
	
	iter++;
      }
    while(goon);
    
    //check that the the routine exited correctly
    check_exit_status_arpack(info);
    
    //free temporary vectors
    nissa_free(temp_x);
    nissa_free(temp_y);
    
    //finalize the eigenstuff calculation
    eigenvalues_of_hermatr_find_arpack_compute(comm,eig_vec,eig_val,neig,min_max,mat_size,target_precision,wspace_size,residue,workd,workl,lworkl,rwork,info,v,ldv);
    
    nissa_free(rwork);
    nissa_free(workl);
    nissa_free(workd);
    nissa_free(v);
    nissa_free(residue);
  }
}

#undef EXTERN_ARPACK
#undef INIT_ARPACK_TO

#endif
