#ifndef _EIGENVALUES_PARPACK_HPP
#define _EIGENVALUES_PARPACK_HPP

#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <mpi.h>

#include <arpack/parpack.hpp>
#include <arpack/debug_c.hpp>

#include "base/vectors.hpp"
#include "eigenvalues/eigenvalues_all.hpp"
#include "new_types/complex.hpp"
#include "routines/ios.hpp"
#include "threads/threads.hpp"

#ifndef EXTERN_PARPACK
 #define EXTERN_PARPACK extern
 #define INIT_PARPACK_TO(...)
#else
 #define INIT_PARPACK_TO(...) __VA_ARGS__
#endif

namespace nissa
{
  EXTERN_PARPACK int use_parpack;
  
  namespace parpack_data
  {
    EXTERN_PARPACK char bmat[2] INIT_PARPACK_TO(="I"); //"I" -> standard eigenvalue problem A*x = lambda*x, "G" -> generalized eigenvalue problem A*x = lambda*M*x
    //largest or smallest
    EXTERN_PARPACK char which[2][3] INIT_PARPACK_TO(={"SM","LM"});
    
    //type for cast
    typedef _Complex double c_t;
    
    class parpack_caller_t
    {
      const int neig;
      const bool min_max;
      const int mat_size;
      const double target_precision;
      const MPI_Fint comm;
      int *glb_info;
      int *glb_ido;
      complex *residue;
      //base for eigenvectors
      const int wspace_size;
      complex *v;
      const int ldv;
      //parameters
      int *iparam;
      //pointer to work
      int *ipntr;
      complex *workd;
      const int lworkl;
      complex *workl;
      complex *rwork;
    public:
      
      //reversed communicator
      int &info=glb_info[0];
      int &ido=glb_ido[0];
      complex *to_be_applied;
      complex *applied;
      
      int nconv() const
      {
	return iparam[4];
      }
      
      parpack_caller_t(const int neig,const bool min_max,const int mat_size,const int mat_size_to_allocate,const double target_precision,const int niter_max,const int wspace_size) :
	neig(neig),
	min_max(min_max),
	mat_size(mat_size),
	target_precision(target_precision),
	comm(MPI_Comm_c2f(MPI_COMM_WORLD)),
	glb_info(nissa_malloc("info",1,int)),
	glb_ido(nissa_malloc("ido",1,int)),
	residue(nissa_malloc("residue",mat_size_to_allocate,complex)),
	wspace_size(wspace_size),
	v(nissa_malloc("v",wspace_size*mat_size_to_allocate,complex)),
	ldv(mat_size_to_allocate),
	iparam(nissa_malloc("iparam",11,int)),
	ipntr(nissa_malloc("ipntr",14,int)),
	workd(nissa_malloc("workd",5*mat_size,complex)),
	lworkl((3*wspace_size+5)*wspace_size),
	workl(nissa_malloc("workl",lworkl,complex)),
	rwork(nissa_malloc("rwork",wspace_size,complex)),
	to_be_applied(nissa_malloc("to_be_applied",mat_size_to_allocate,complex)),
	applied(nissa_malloc("applied",mat_size_to_allocate,complex))
      {
	info=0;
	ido=0;

#ifdef ENABLE_PARPACK_DEBUG
	// see debug.doc in the arpack documentation for the mcXYZ parameters
	const int logfil=6;   // standard output
	const int ndigit=-10; // 3 digits and 72 columns (positive would print on 132 columns)
	const int mgetv0=0;   // do not print residual vector
	const int mcaupd=1;
	const int mcaupd2=1;
	const int mcaitr=1;
	const int mceigh=0;
	const int mcgets=0;
	const int mcapps=0;
	const int mceupd=1;
	
	debug_c(logfil,ndigit,mgetv0,
          0,0,0,0,0,0,0,
          0,0,0,0,0,0,0,
          mcaupd,mcaupd2,mcaitr,mceigh,mcgets,mcapps,mceupd);
#endif
	
	//init params
	vector_reset(ipntr);
	vector_reset(iparam);
	iparam[0]=1;
	iparam[2]=niter_max;
	iparam[3]=1;
	iparam[6]=1;
	THREAD_BARRIER();
	
	//print some debug info
	VERBOSITY_LV2_MASTER_PRINTF("mat_size: %d\n",mat_size);
	VERBOSITY_LV2_MASTER_PRINTF("which[min_max]: %s\n",which[min_max]);
	VERBOSITY_LV2_MASTER_PRINTF("neig: %d\n",neig);
	VERBOSITY_LV2_MASTER_PRINTF("target_precision: %lg\n",target_precision);
	VERBOSITY_LV2_MASTER_PRINTF("wspace_size: %d\n",wspace_size);
	VERBOSITY_LV2_MASTER_PRINTF("ldv: %d\n",ldv);
	for(int i=0;i<11;i++)
	  VERBOSITY_LV2_MASTER_PRINTF("iparam[%d]: %d\n",i,iparam[i]);
	for(int i=0;i<14;i++)
	  VERBOSITY_LV2_MASTER_PRINTF("ipntr[%d]: %d\n",i,ipntr[i]);
	VERBOSITY_LV2_MASTER_PRINTF("lworkl: %d\n",lworkl);
	VERBOSITY_LV2_MASTER_PRINTF("info: %d\n",info);
      }
      
#define REST_OF_PARS		\
      bmat,			\
	mat_size,		\
	which[min_max],		\
	neig,			\
	target_precision,	\
	(c_t*)residue,		\
	wspace_size,		\
	(c_t*)v,		\
	ldv,			\
	iparam,			\
	ipntr,			\
	(c_t*)workd,		\
	(c_t*)workl,		\
	lworkl,			\
	(c_t*)rwork,		\
	&info
      
      int iteration()
      {
	
	//If not first iteration, read the result from the reversed communicator
	if(ido!=0)
	  {
	    complex *y=workd+ipntr[1]-1;
	    memcpy(y,applied,sizeof(complex)*mat_size);
	  }
	
	//Only master thread
	THREAD_BARRIER();
	if(IS_MASTER_THREAD)
	  arpack::internal::pznaupd_c(comm,             //MPI handle
				      &ido,             //Reverse communication flag
				      REST_OF_PARS);
	THREAD_BARRIER();
	
	//decrypt info
	switch(info)
	  {
	  case 0: /*Normal execution */;break;
	  case 1: CRASH("Maximum number of iterations taken.");break;
	  case 3: CRASH("No shifts could be applied during a cycle of the "
			"Implicitly restarted Arnoldi iteration. One possibility "
			"is to increase the size of NCV relative to NEV.");break;
	  case -1: CRASH("N must be positive.");break;
	  case -2: CRASH("NEV must be positive.");break;
	  case -3: CRASH("NCV-NEV >= 2 and less than or equal to N.");break;
	  case -4: CRASH("The maximum number of Arnoldi update iteration must be greater than zero.");break;
	  case -5: CRASH("WHICH must be one of 'LM', 'SM', 'LR', 'SR', 'LI', 'SI'");break;
	  case -6: CRASH("BMAT must be one of 'I' or 'G'.");break;
	  case -7: CRASH("Length of private work array is not sufficient.");break;
	  case -8: CRASH("Error return from LAPACK eigenvalue calculation;");break;
	  case -9: CRASH("Starting vector is zero.");break;
	  case -10: CRASH("IPARAM(7) must be 1,2,3.");break;
	  case -11: CRASH("IPARAM(7) = 1 and BMAT = 'G' are incompatable.");break;
	  case -12: CRASH("IPARAM(1) must be equal to 0 or 1.");break;
	  case -9999: CRASH("Could not build an Arnoldi factorization.");break;
	  default:crash("Unknown error code");break;
	  }
	
	//perform the required task
	int goon;
	switch(ido)
	  {
	    complex *x;
	  case 1:
	    goon=true;
	    
	    //reverse communication
	    x=workd+ipntr[0]-1;
	    memcpy(to_be_applied,x,sizeof(complex)*mat_size);
	    set_borders_invalid(to_be_applied);
	    break;
	  case 99:
	    goon=false;
	    break;
	  default:
	    goon=false;
	    CRASH("ido %d not contemplated",ido);
	  }
	
	return goon;
      }
      
      void finalize(complex **eig_vec,complex *eig_val)
      {
	
	THREAD_BARRIER();
	
	if(IS_MASTER_THREAD)
	  {
	    //parameters for eigenvector calculation
	    complex *ceig_val=new complex[neig+1];
	    int *select=new int[wspace_size];
	    complex *workev=new complex[2*wspace_size];
	    complex *temp_vec=new complex[neig*mat_size];
	    
	    arpack::internal::pzneupd_c(comm,             //Reverse communicator
					1,                //1: Compute Ritz vectors or Schur vectors. 0: Ritz values only
					"A",		//'A': Compute NEV Ritz vectors; 'P': Compute NEV Schur vectors; 'S': compute some of the Ritz vectors
					select,		//Logical array of dimension wspace_size
					(c_t*)ceig_val,	//Complex array of dimension neig+1.
					(c_t*)temp_vec,	//Complex mat_size by neig array
					mat_size,		//The leading dimension of the array temp_vec
					0.0,		//If IPARAM(7) = 3 then represents the shift
					(c_t*)workev,     //Complex*16 work array of dimension 2*NCV
					REST_OF_PARS);
	    
#undef REST_OF_PARS
	    
	    switch(info)
	      {
	      case 0: /* Normal exit */;break;
	      case 1: CRASH("The Schur form computed by LAPACK routine csheqr could not be reordered by LAPACK routine ztrsen. "
			    "Re-enter subroutine zneupd with IPARAM(5) = NCV and increase the size of the array D to have dimension at least dimension NCV and allocate at least NCV columns for Z."
			    "NOTE: Not necessary if Z and V share the same space. Please notify the authors if this error occurs.");break;
	      case -1: CRASH("N=mat_size=%d must be positive.",mat_size);break;
	      case -2: CRASH("NEV=neig=%d must be positive.",neig);break;
	      case -3: CRASH("NCV-NEV=wspace_size-neig=%d-%d=%d >= 1 and less than or equal to N=mat_size=%d.",wspace_size,neig,wspace_size-neig,mat_size);break;
	      case -5: CRASH("WHICH must be one of 'LM', 'SM', 'LR', 'SR', 'LI', 'SI', is '%s'.",which[min_max]);break;
	      case -6: CRASH("BMAT must be one of 'I' or 'G', is '%s'.",bmat);break;
	      case -7: CRASH("Length of private work WORKL array=%d is not sufficient.",lworkl);break;
	      case -8: CRASH("Error return from LAPACK eigenvalue calculation. This should never happened.");break;
	      case -9: CRASH("Error return from calculation of eigenvectors. Informational error from LAPACK routine ztrevc.");break;
	      case -10: CRASH("IPARAM(7) must be 1, 2, 3.");break;
	      case -11: CRASH("IPARAM(7) = 1 and BMAT = 'G' are incompatible.");break;
	      case -12: CRASH("HOWMANY = 'S' not yet implemented.");break;
	      case -13: CRASH("HOWMANY must be one of 'A' or 'P' if RVEC = .true.");break;
	      case -14: CRASH("ZNAUPD did not find any eigenvalues to sufficient accuracy.");break;
	      case -15: CRASH("ZNEUPD got a different count of the number of converged Ritz values than ZNAUPD got."
			      "This indicates the user probably made an error in passing data from ZNAUPD to ZNEUPD or that the data was modified before entering ZNEUPD.");break;
	      }
	    
	    //find the ordering
	    std::vector<std::pair<double,int>> ord;
	    for(int i=0;i<neig;i++)
	      ord.push_back({complex_norm2(ceig_val[i]),i});
	    std::sort(ord.begin(),ord.end());
	    
	    //store result
	    for(int i=0;i<neig;i++)
	      {
		int iin=ord[i].second;
		memcpy(eig_vec[i],temp_vec+mat_size*iin,sizeof(complex)*mat_size);
		complex_copy(eig_val[i],ceig_val[iin]);
	      }
	    
	    delete[] temp_vec;
	    delete[] workev;
	    delete[] ceig_val;
	    delete[] select;
	  }
	
	THREAD_BARRIER();
	
	//Ensures that all threads have the same eig_vec
	complex *master_eig_val;
	THREAD_BROADCAST_PTR(master_eig_val,eig_val);
	if(master_eig_val!=eig_val)
	  for(int i=0;i<neig;i++)
	    complex_copy(eig_val[i],master_eig_val[i]);
	THREAD_BARRIER();
      }
      
      ~parpack_caller_t()
      {
	nissa_free(glb_info);
	nissa_free(glb_ido);
	nissa_free(residue);
	nissa_free(v);
	nissa_free(iparam);
	nissa_free(ipntr);
	nissa_free(workd);
	nissa_free(workl);
	nissa_free(rwork);
	nissa_free(to_be_applied);
	nissa_free(applied);
      }
    };
  }
  
  //use parpack
  template <class Fmat,class Filler>
  void eigenvalues_find_parpack(complex **eig_vec,complex *eig_val,int neig,bool min_max,
				const int mat_size,const int mat_size_to_allocate,const Fmat &imp_mat,
				const double target_precision,const int niter_max,
				const Filler &filler,int wspace_size=DEFAULT_EIGPROB_WSPACE_SIZE)
  {
    using namespace parpack_data;
    wspace_size=std::min(std::max(2*neig,wspace_size),mat_size);
    parpack_caller_t caller(neig,min_max,mat_size,mat_size_to_allocate,target_precision,niter_max,wspace_size);
    
    //main loop to find Ritz basis
    bool goon;
    int iter=0;
    do
      {
	//invoke the step
	goon=caller.iteration();
#ifndef ENABLE_PARPACK_DEBUG
	VERBOSITY_LV1_MASTER_PRINTF("iteration %d, ido: %d, info: %d\n",iter,caller.ido,caller.info);
#endif	
	if(goon)
	  imp_mat(caller.applied,caller.to_be_applied);
	
	iter++;
      }
    while(goon and iter<niter_max);
    
    caller.finalize(eig_vec,eig_val);
  }
}

#undef EXTERN_PARPACK
#undef INIT_PARPACK_TO

#endif
