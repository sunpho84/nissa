#ifdef HAVE_CONFIG_H
# include "config.hpp"
#endif

#define EXTERN_DEBUG
# include "base/debug.hpp"

#include <signal.h>
#include <execinfo.h>
#ifdef USE_MPI
# include <mpi.h>
#endif
#include <stdarg.h>
#include <stdlib.h>
#include <unistd.h>

#ifdef USE_CUDA
# include <cuda_runtime.h>
#endif

#include "base/field.hpp"
#include "geometry/geometry_lx.hpp"
#include "new_types/float_128.hpp"
#include "routines/ios.hpp"
#include "routines/mpi_routines.hpp"
#include "threads/threads.hpp"

#include "vectors.hpp"

namespace nissa
{
  /// Implements the trap to debug
  void debug_loop()
  {
    if(is_master_rank())
      {
	volatile int flag=0;
	
	printf("Entering debug loop on rank %d, flag has address %p please type:\n"
	       "$ gdb -p %d\n"
	       "$ set flag=1\n"
	       "$ continue\n",
	       rank,
	       &flag,
	       getpid());
	
	while(flag==0);
      }
    
    ranks_barrier();
  }
  
  //take the time
  double take_time()
  {
#ifdef USE_MPI
    return MPI_Wtime();
#else
    return (double)clock()/CLOCKS_PER_SEC;
#endif
  }
  
  //write the list of called routines
  void print_backtrace_list(int which_rank)
  {
    void *callstack[128];
    int frames=backtrace(callstack,128);
    char **strs=backtrace_symbols(callstack,frames);
    
    if(which_rank==-1 or which_rank==rank)
      {
	printf("Backtracing...\n");
	for(int i=0;i<frames;i++) printf("%s\n",strs[i]);
      }
    
    free(strs);
  }
  
  //crash reporting the expanded error message
  CUDA_HOST_AND_DEVICE
  void internal_crash(int line,const char *file,const char *templ,...)
  {
#ifndef COMPILING_FOR_DEVICE
    fflush(stdout);
    fflush(stderr);
    
    //give time to master thread to crash, if possible
    sleep(1);
    
    // if(is_master_rank())
      {
	//expand error message
	char mess[1024];
	va_list ap;
	va_start(ap,templ);
	vsprintf(mess,templ,ap);
	va_end(ap);
	
	fprintf(stderr,RED_HIGHLIGHT "on rank %d ERROR on line %d of file \"%s\", message error: \"%s\"." DO_NOT_HIGHLIGHT,rank,line,file,mess);
	fprintf(stderr,"Memory used: %ld bytes per rank (%ld bytes total)\n",required_memory,required_memory*nranks);
	print_backtrace_list();
    }
    
    ranks_abort(0);
#else
    __trap();
#endif
  }
  
  void internal_crash_printing_error(int line,const char *file,int err_code,const char *templ,...)
  {
    if(err_code)
      {
	//print error code
	char str1[1024];
	sprintf(str1,"returned code %d",err_code);
	
	//expand error message
	char str2[1024];
	va_list ap;
	va_start(ap,templ);
	vsprintf(str2,templ,ap);
	va_end(ap);
	
	internal_crash(line,file,"%s %s",str1,str2);
      }
  }
  
  //called when signal received
  void signal_handler(int sig)
  {
    MASTER_PRINTF("maximal memory used: %ld\n",max_required_memory);
    verbosity_lv=3;
    char name[100];
    switch(sig)
      {
      case SIGSEGV: sprintf(name,"segmentation violation");break;
      case SIGFPE: sprintf(name,"floating-point exception");break;
      case SIGXCPU: sprintf(name,"cpu time limit exceeded");break;
      case SIGBUS: sprintf(name,"bus error");break;
      case SIGINT: sprintf(name," program interrupted");break;
      case SIGABRT: sprintf(name,"abort signal");break;
      default: sprintf(name,"unassociated");break;
      }
    print_backtrace_list();
    print_all_vect_content();
    CRASH("signal %d (%s) detected, exiting",sig,name);
  }
  
#ifdef USE_MPI
  //decript the MPI error
  void internal_decrypt_MPI_error(int line,const char *file,int rc,const char *templ,...)
  {
    if(rc!=MPI_SUCCESS and is_master_rank())
      {
	char err[1024];
	int len=1024;
	MPI_Error_string(rc,err,&len);
	
	va_list ap;
	va_start(ap,templ);
	char mess[1024];
	vsprintf(mess,templ,ap);
	va_end(ap);
	
	internal_crash(line,file,"%s, MPI raised error: %s",mess,err);
      }
  }
#endif
  
#if USE_CUDA
  void internal_decrypt_cuda_error(int line,const char *file,cudaError_t rc,const char *templ,...)
  {
    if(rc!=cudaSuccess and rank==0)
      {
	va_list ap;
	va_start(ap,templ);
	char mess[1024];
	vsprintf(mess,templ,ap);
	va_end(ap);
	
	internal_crash(line,file,"%s, cuda raised error: %s",mess,cudaGetErrorString(rc));
      }
  }
#endif
  
  //perform a simple check on 128 bit precision
  void check_128_bit_prec()
  {
    Float128 a=1;
    a+=1e-20;
    a+=-1;
    
    const double res=a.roundDown();
    if(fabs(res-1e-20)>1e-30)
      CRASH("float_128, 1+1e-20-1=%lg, difference with 1e-20: %lg",res,res-1e-20);
    VERBOSITY_LV2_MASTER_PRINTF("128 bit precision is working, 1+1e-20-1=%lg where %lg expected in double prec\n",res,1+1e-20-1);
  }
  
  std::string siteAsString(const int& n)
  {
    const auto c=glbCoordOfGlblx(n);
    
    std::string res=
      std::to_string(n)+
      "("+std::to_string(c[0]);
    
    for(int nu=1;nu<NDIM;nu++)
      res+=","+std::to_string(c[nu]);
    res+=")";
    
    return res;
  }
  
  /// Spoil the receiving and sending buffer to check consistency
  void taintTheCommBuffers()
  {
    PAR(0,recv_buf_size/sizeof(int),
	CAPTURE(),
	i,
	{
	  ((int*)recv_buf)[i]=-3;
	});
    
    PAR(0,send_buf_size/sizeof(int),
	CAPTURE(),
	i,
	{
	  ((int*)send_buf)[i]=-4;
	});
  }
  
  void testLxHaloExchange()
  {
    LxField<int> test("testHalo",WITH_HALO);
    
    PAR(0,locVol,
	CAPTURE(TO_WRITE(test)),i,
	{
	  test[i]=glblxOfLoclx[i];
	});
    
    test.invalidateHalo();
    
    taintTheCommBuffers();
    
    test.updateHalo();
    
    NISSA_LOC_VOL_LOOP(site)
      {
	for(int ori=0;ori<2;ori++)
	  for(int mu=0;mu<NDIM;mu++)
	    {
	      const int ln=loclx_neigh[ori][site][mu];
	      const int gn=(ln<locVol)?glblxOfLoclx[ln]:glblxOfBordlx[ln-locVol];
	      const int neighVal=test[ln];
	      
	      if(neighVal!=gn)
		MASTER_PRINTF("site %s ori %d dir %d has neigh %s with val %s\n",
			      siteAsString(glblxOfLoclx[site]).c_str(),
			      ori,mu,
			      siteAsString(gn).c_str(),
			      siteAsString(neighVal).c_str());
	    }
      }
    
    MASTER_PRINTF("lx halo communicates consistently\n");
  }
  
  void testEoHaloExchange()
  {
    EoField<int> test("testEoHalo",WITH_HALO);
    
    FOR_BOTH_PARITIES(par,
		      PAR(0,locVolh,
			  CAPTURE(par,test=test[par].getWritable()),i,
			  {
			    test[i]=glblxOfLoclx[loclx_of_loceo[par][i]];
			  });
		      );
    
    test.invalidateHalo();
    
    taintTheCommBuffers();
    
    test.updateHalo();
    
    FOR_BOTH_PARITIES(par,
		      for(int eoSite=0;eoSite<locVolh;eoSite++)
		      {
			for(int ori=0;ori<2;ori++)
			  for(int mu=0;mu<NDIM;mu++)
			    {
			      const int ln=((ori==0)?loceo_neighdw:loceo_neighup)[par][eoSite][mu];
			      const int lx=loclx_of_loceo[!par][ln];
			      const int gn=(lx<locVol)?glblxOfLoclx[lx]:glblxOfBordlx[lx-locVol];
			      const int neighVal=test[!par][ln];
			      
			      if(neighVal!=gn)
				MASTER_PRINTF("site %s ori %d dir %d has neigh %s with val %s\n",
					      siteAsString(glblxOfLoclx[loclx_of_loceo[!par][eoSite]]).c_str(),
					      ori,mu,
					      siteAsString(gn).c_str(),
					      siteAsString(neighVal).c_str());
			    }
		      });
  
    MASTER_PRINTF("eo edges communicates consistently\n");
  }
  
  void testLxEdgesExchange()
  {
    LxField<int> test("testEdge",WITH_HALO_EDGES);
    PAR(0,locVol,
	CAPTURE(TO_WRITE(test)),i,
	{
	  test[i]=glblxOfLoclx[i];
	});
    
    PAR(0,bordVol,
	CAPTURE(TO_WRITE(test)),i,
	{
	  test[i+locVol]=-1;
	});
    
    PAR(0,edgeVol,
	CAPTURE(TO_WRITE(test)),i,
	{
	  test[i+locVol+bordVol]=-2;
	});
    
    test.invalidateHalo();
    test.invalidateEdges();
    test.updateHalo();
    
    taintTheCommBuffers();
    
    test.updateEdges();
    
    NISSA_LOC_VOL_LOOP(site)
      {
	for(int ori1=0;ori1<2;ori1++)
	  for(int ori2=0;ori2<2;ori2++)
	    for(int iEdge=0;iEdge<nEdges;iEdge++)
	      {
		const auto [mu,nu]=edge_dirs[iEdge];
		
		const int l1n=loclx_neigh[ori1][site][mu];
		const int ln=loclx_neigh[ori2][l1n][nu];
		const int gn=(ln<locVol)?glblxOfLoclx[ln]:((ln<locVol+bordVol)?glblxOfBordlx[ln-locVol]:glblxOfEdgelx[ln-locVol-bordVol]);
		const int neighVal=test[ln];
		
		if(neighVal!=gn)
		  MASTER_PRINTF("site %s ori (%d,%d) dir (%d,%d) has edgelx neigh %s with val %s\n",
				siteAsString(glblxOfLoclx[site]).c_str(),
				ori1,ori2,mu,nu,
				siteAsString(gn).c_str(),
				siteAsString(neighVal).c_str());
	      }
      }
    MASTER_PRINTF("lx edges communicates consistently\n");
  }
  
  void testEoEdgesExchange()
  {
    EoField<int> test("testEoEdge",WITH_HALO_EDGES);
    FOR_BOTH_PARITIES(par,
		      PAR(0,locVolh,
			  CAPTURE(par,test=test[par].getWritable()),i,
			  {
			    test[i]=glblxOfLoclx[loclx_of_loceo[par][i]];
			  });
		      
		      PAR(0,bordVolh,
			  CAPTURE(test=test[par].getWritable()),i,
			  {
			    test[i+locVolh]=-1;
			  });
		      
		      PAR(0,edgeVolh,
			  CAPTURE(test=test[par].getWritable()),i,
			  {
			    test[i+locVolh+bordVolh]=-2;
			  });
		      );
    
    int* r=(int*)recv_buf;
    int* s=(int*)send_buf;
    const int n=recv_buf_size/sizeof(int);
    PAR(0,n,
	CAPTURE(r,s),
	i,
	{
	  r[i]=-3;
	  s[i]=-4;
	});
    
    test.invalidateHalo();
    test.invalidateEdges();
    test.updateHalo();
    test.updateEdges();
    
    FOR_BOTH_PARITIES(par,
		      NISSA_LOC_VOLH_LOOP(ieo)
		      {
			for(int ori1=0;ori1<2;ori1++)
			  for(int ori2=0;ori2<2;ori2++)
			    for(int iEdge=0;iEdge<nEdges;iEdge++)
			      {
				const auto [mu,nu]=edge_dirs[iEdge];
				
				const int l1n=((ori1==0)?loceo_neighdw:loceo_neighup)[par][ieo][mu];
				const int ln=((ori2==0)?loceo_neighdw:loceo_neighup)[!par][l1n][nu];
				const int lx=loclx_of_loceo[par][ln];
				const int gn=(lx<locVol)?glblxOfLoclx[lx]:((lx<locVol+bordVol)?glblxOfBordlx[lx-locVol]:glblxOfEdgelx[lx-locVol-bordVol]);
				const int neighVal=test[par][ln];
				
				if(neighVal!=gn)
				  MASTER_PRINTF("par %d site %s ori (%d,%d) dir (%d,%d) neigh %s with val %s\n",
						par(),
						siteAsString(glblxOfLoclx[loclx_of_loceo[par][ieo]]).c_str(),
						ori1,ori2,mu,nu,
						siteAsString(gn).c_str(),
						siteAsString(neighVal).c_str());
			      }
		      }
		      );
    
    MASTER_PRINTF("eo edges communicates consistently\n");
  }
}
