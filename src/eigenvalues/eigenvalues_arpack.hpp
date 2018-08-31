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
    const char which[2][3]={"LM","SM"};
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
    
    //type for cast
    typedef _Complex double c_t;
    
    //main loop to find Ritz basis
    int info=0;
    int ido=0;
    do
      {
	//invoke the step
	arpack::internal::pznaupd_c(comm,&ido,bmat,mat_size,which[min_max],neig,target_precision,(c_t*)residue,wspace_size,(c_t*)v,ldv,iparam,ipntr,(c_t*)workd,(c_t*)workl,lworkl,(c_t*)rwork,&info);
	verbosity_lv1_master_printf("ido: %d, info: %d, ipntr[0]: %d, ipntr[1]: %d\n",ido,info,ipntr[0],ipntr[1]);
	
	//reverse communication
	complex *x=workd+ipntr[0]-1;
	complex *y=workd+ipntr[1]-1;
	imp_mat(y,x);
      }
    while(ido==1);
    
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
