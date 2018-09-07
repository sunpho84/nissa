#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#define EXTERN_ARPACK
#include "eigenvalues_arpack.hpp"

namespace nissa
{
  //flag to use arpack or not
  int use_arpack;
  
  namespace arpack_data
  {
    //check the exit status
    void check_exit_status_arpack(const int info)
    {
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
    }
    
    //initialize parameters
    void init_params(int *iparam,int niter_max)
    {
      //init params
      iparam[0]=1;
      iparam[2]=niter_max;
      iparam[3]=1;
      iparam[6]=1;
    }
    
    //common chunk to compute eigenvalues on the basis of Ritz vectors
    void eigenvalues_of_hermatr_find_arpack_compute(MPI_Fint &comm,complex **eig_vec,complex *eig_val,int neig,bool min_max,
						    const int mat_size,const double target_precision,const int wspace_size,
						    complex *residue,complex *workd,complex *workl,const int lworkl,complex *rwork,int &info,
						    complex *v,const int ldv)
    {
      using namespace arpack_data;
      
      //parameters for eigenvector calculation
      const int rvec=1;
      const char howmny[]="P";
      complex ceig_val[neig];
      const c_t sigma=0.0;
      int select[wspace_size];
      complex *workev=nissa_malloc("workev",2*wspace_size,complex);
      complex *temp_vec=nissa_malloc("temp_vec",neig*mat_size,complex);
      arpack::internal::pzneupd_c(comm,rvec,howmny,select,(c_t*)ceig_val,(c_t*)temp_vec,mat_size,sigma,(c_t*)workev,bmat,mat_size,which[min_max],neig,target_precision,(c_t*)residue,wspace_size,(c_t*)v,ldv,iparam,ipntr,(c_t*)workd,(c_t*)workl,lworkl,(c_t*)rwork,&info);
      
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
      
      nissa_free(temp_vec);
      nissa_free(workev);
    }
  }
}
