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

namespace nissa
{
  //type for cast
  typedef _Complex double c_t;
  
  //use arpack
  template <class Fmat,class Filler>
  void eigenvalues_of_hermatr_find_arpack(complex **eig_vec,double *eig_val,int neig,bool min_max,
					  const int mat_size,const int mat_size_to_allocate,const Fmat &imp_mat,
					  const double target_precision,const int niter_max,
					  const Filler &filler)
  {
    //convert the communicator
    MPI_Fint comm=MPI_Comm_c2f(MPI_COMM_WORLD);
    
    //normal eigenvalue problem
    const char bmat[]="I";
    //largest or smallest
    const char which[2][3]={"SM","LM"};
    //vector of residue
    complex *residue=nissa_malloc("residue",mat_size_to_allocate,complex);
    //base for eigenvectors
    const int wspace_size=2*neig;
    //workspace
    complex *v=nissa_malloc("v",wspace_size*mat_size_to_allocate,complex);
    //size of each line
    const int ldv=mat_size_to_allocate;
    //parameters
    int iparam[11];
    //pointer to work
    int ipntr[14];
    //internal workspace
    complex *workd=nissa_malloc("workd",3*mat_size,complex);
    int lworkl=(3*wspace_size+5)*wspace_size;
    complex *workl=nissa_malloc("workl",lworkl,complex);
    complex *rwork=nissa_malloc("rwork",wspace_size,complex);
    
    //init params
    iparam[0]=1;
    iparam[2]=niter_max;
    iparam[3]=1;
    iparam[6]=1;
    
    //temporary vectors
    complex *temp_x=nissa_malloc("temp_x",mat_size_to_allocate,complex);
    complex *temp_y=nissa_malloc("temp_y",mat_size_to_allocate,complex);
    
    //main loop to find Ritz basis
    int info=0;
    int ido=0;
    int iter=0;
    do
      {
	//invoke the step
	arpack::internal::pznaupd_c(comm,&ido,bmat,mat_size,which[min_max],neig,target_precision,(c_t*)residue,wspace_size,(c_t*)v,ldv,iparam,ipntr,(c_t*)workd,(c_t*)workl,lworkl,(c_t*)rwork,&info);
	verbosity_lv1_master_printf("iteration %d, ido: %d, info: %d, nconv: %d/%d\n",iter,ido,info,ipntr[0],ipntr[1],iparam[4],neig);
	
	//reverse communication
	complex *x=workd+ipntr[0]-1;
	complex *y=workd+ipntr[1]-1;
	
	memcpy(temp_x,x,sizeof(complex)*mat_size);
	set_borders_invalid(temp_x);
	imp_mat(temp_y,temp_x);
	memcpy(y,temp_y,sizeof(complex)*mat_size);
	
	iter++;
      }
    while(ido==1);
    
    //decrypt info
    switch(info)
      {
      case 0: master_printf("Normal exit.\n");break;
      case 1: crash("Maximum number of iterations taken.");break;
      case 3: crash("No shifts could be applied during a cycle of the "
			    "Implicitly restarted Arnoldi iteration. One possibility "
			    "is to increase the size of NCV relative to NEV.");break;
      case-1: crash("N must be positive.");break;
      case-2: crash("NEV must be positive.");break;
      case-3: crash("NCV-NEV >= 2 and less than or equal to N.");break;
      case-4: crash("The maximum number of Arnoldi update iteration must be greater than zero.");break;
      case-5: crash("WHICH must be one of 'LM', 'SM', 'LR', 'SR', 'LI', 'SI'");break;
      case-6: crash("BMAT must be one of 'I' or 'G'.");break;
      case-7: crash("Length of private work array is not sufficient.");break;
      case-8: crash("Error return from LAPACK eigenvalue calculation;");break;
      case-9: crash("Starting vector is zero.");break;
      case-10: crash("IPARAM(7) must be 1,2,3.");break;
      case-11: crash("IPARAM(7) = 1 and BMAT = 'G' are incompatable.");break;
      case-12: crash("IPARAM(1) must be equal to 0 or 1.");break;
      case-9999: crash("Could not build an Arnoldi factorization.");break;
      default:crash("Unknown error code");break;
      }
    
    //free temporary vectors
    nissa_free(temp_x);
    nissa_free(temp_y);
    
    //parameters for eigenvector calculation
    const int rvec=1;
    const char howmny[]="P";
    complex ceig_val[neig];
    const c_t sigma=0.0;
    int select[wspace_size];
    complex *workev=nissa_malloc("workev",2*wspace_size,complex);
    complex *temp_vec=nissa_malloc("temp_vec",neig*mat_size,complex);
    arpack::internal::pzneupd_c(comm,rvec,howmny,select,(c_t*)ceig_val,(c_t*)temp_vec,mat_size,sigma,(c_t*)workev,bmat,mat_size,which[min_max],neig,target_precision,(c_t*)residue,wspace_size,(c_t*)v,ldv,iparam,ipntr,(c_t*)workd,(c_t*)workl,lworkl,(c_t*)rwork,&info);
    
    //store result
    for(int i=0;i<neig;i++)
      {
	memcpy(eig_vec[i],temp_vec+mat_size*i,sizeof(complex)*mat_size);
	eig_val[i]=ceig_val[i][RE];
      }
    
    
    nissa_free(temp_vec);
    nissa_free(workev);
    nissa_free(rwork);
    nissa_free(workl);
    nissa_free(workd);
    nissa_free(v);
    nissa_free(residue);
  }
} 

#endif
