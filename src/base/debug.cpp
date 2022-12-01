#ifdef HAVE_CONFIG_H
# include "config.hpp"
#endif

#define EXTERN_DEBUG
# include "base/debug.hpp"

#include <signal.h>
#include <errno.h>
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
  void print_backtrace_list()
  {
    void *callstack[128];
    int frames=backtrace(callstack,128);
    char **strs=backtrace_symbols(callstack,frames);
    
    //only master rank, but not master thread
    if(is_master_rank())
      {
	printf("Backtracing...\n");
	for(int i=0;i<frames;i++) printf("%s\n",strs[i]);
      }
    
    free(strs);
  }
  
  //crash reporting the expanded error message
  void internal_crash(int line,const char *file,const char *templ,...)
  {
    fflush(stdout);
    fflush(stderr);
    
    //give time to master thread to crash, if possible
    if(not IS_MASTER_THREAD)
      sleep(1);
    
    if(is_master_rank())
      {
	//expand error message
	char mess[1024];
	va_list ap;
	va_start(ap,templ);
	vsprintf(mess,templ,ap);
	va_end(ap);
	
	fprintf(stderr,"\x1b[31m" "ERROR on line %d of file \"%s\", message error: \"%s\".\n\x1b[0m",line,file,mess);
	fprintf(stderr,"Memory used: %ld bytes per rank (%ld bytes total)\n",required_memory,required_memory*nranks);
	print_backtrace_list();
	ranks_abort(0);
      }
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
    master_printf("maximal memory used: %ld\n",max_required_memory);
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
    crash("signal %d (%s) detected, exiting",sig,name);
  }
  
#ifdef USE_MPI
  //decript the MPI error
  void internal_decript_MPI_error(int line,const char *file,int rc,const char *templ,...)
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
  void internal_decript_cuda_error(int line,const char *file,cudaError_t rc,const char *templ,...)
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
    float_128 a;
    float_128_from_64(a,1);
    float_128_summassign_64(a,1e-20);
    float_128_summassign_64(a,-1);
    
    double res=a[0]+a[1];
    if(fabs(res-1e-20)>1e-30) crash("float_128, 1+1e-20-1=%lg, difference with 1e-20: %lg",res,res-1e-20);
    verbosity_lv2_master_printf("128 bit precision is working, 1+1e-20-1=%lg where %lg expected in double prec\n",res,1+1e-20-1);
  }
  
  std::string siteAsString(const int& n)
  {
    const auto c=glb_coord_of_glblx(n);
    
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
    int* r=(int*)recv_buf;
    NISSA_PARALLEL_LOOP(i,0,recv_buf_size/sizeof(int))
      r[i]=-3;
    NISSA_PARALLEL_LOOP_END;
    
    int* s=(int*)send_buf;
    NISSA_PARALLEL_LOOP(i,0,send_buf_size/sizeof(int))
      s[i]=-4;
    NISSA_PARALLEL_LOOP_END;
  }
  
  void testLxHaloExchange()
  {
    LxField<int> test("testHalo",WITH_HALO);
    NISSA_PARALLEL_LOOP(i,0,locVol)
      test[i]=glblxOfLoclx[i];
    NISSA_PARALLEL_LOOP_END;
    
    test.invalidateHalo();
    
    taintTheCommBuffers();
    
    test.updateHalo();
    
    NISSA_PARALLEL_LOOP(site,0,locVol)
      {
	for(int ori=0;ori<2;ori++)
	  for(int mu=0;mu<NDIM;mu++)
	    {
	      const int ln=loclx_neigh[ori][site][mu];
	      const int gn=(ln<locVol)?glblxOfLoclx[ln]:glblxOfBordlx[ln-locVol];
	      const int neighVal=test[ln];
	      
	      if(neighVal!=gn)
		master_printf("site %s ori %d dir %d has neigh %s with val %s\n",
			      siteAsString(glblxOfLoclx[site]).c_str(),
			      ori,mu,
			      siteAsString(gn).c_str(),
			      siteAsString(neighVal).c_str());
	    }
      }
    NISSA_PARALLEL_LOOP_END;
    
    master_printf("lx halo communicates consistently\n");
  }
  
  void testEoHaloExchange()
  {
    EoField<int> test("testEoHalo",WITH_HALO);
    forBothParities([&test](const auto& par)
    {
      NISSA_PARALLEL_LOOP(i,0,locVolh)
	test[par][i]=glblxOfLoclx[loclx_of_loceo[par][i]];
      NISSA_PARALLEL_LOOP_END;
    });
    test.invalidateHalo();
    
    taintTheCommBuffers();
    
    test.updateHalo();
    
    forBothParities([&test](const auto& par)
    {
      NISSA_PARALLEL_LOOP(eoSite,0,locVolh)
	{
	  for(int ori=0;ori<2;ori++)
	    for(int mu=0;mu<NDIM;mu++)
	      {
		const int ln=((ori==0)?loceo_neighdw:loceo_neighup)[par][eoSite][mu];
		const int lx=loclx_of_loceo[!par][ln];
		const int gn=(lx<locVol)?glblxOfLoclx[lx]:glblxOfBordlx[lx-locVol];
		const int neighVal=test[!par][ln];
		
		if(neighVal!=gn)
		  master_printf("site %s ori %d dir %d has neigh %s with val %s\n",
				siteAsString(glblxOfLoclx[loclx_of_loceo[!par][eoSite]]).c_str(),
				ori,mu,
				siteAsString(gn).c_str(),
				siteAsString(neighVal).c_str());
	      }
	}
      NISSA_PARALLEL_LOOP_END;
    });
    
    master_printf("eo edges communicates consistently\n");
  }
  
  void testLxEdgesExchange()
  {
    LxField<int> test("testEdge",WITH_HALO_EDGES);
    NISSA_PARALLEL_LOOP(i,0,locVol)
      test[i]=glblxOfLoclx[i];
    NISSA_PARALLEL_LOOP_END;
    
    NISSA_PARALLEL_LOOP(i,0,bord_vol)
      test[i+locVol]=-1;
    NISSA_PARALLEL_LOOP_END;
    
    NISSA_PARALLEL_LOOP(i,0,edge_vol)
      test[i+locVol+bord_vol]=-2;
    NISSA_PARALLEL_LOOP_END;
    
    test.invalidateHalo();
    test.invalidateEdges();
    test.updateHalo();
    
    taintTheCommBuffers();
    
    test.updateEdges();
    
    NISSA_PARALLEL_LOOP(site,0,locVol)
      {
	for(int ori1=0;ori1<2;ori1++)
	  for(int ori2=0;ori2<2;ori2++)
	    for(int iEdge=0;iEdge<nEdges;iEdge++)
	      {
		const auto [mu,nu]=edge_dirs[iEdge];
		
		const int l1n=loclx_neigh[ori1][site][mu];
		const int ln=loclx_neigh[ori2][l1n][nu];
		const int gn=(ln<locVol)?glblxOfLoclx[ln]:((ln<locVol+bord_vol)?glblxOfBordlx[ln-locVol]:glblxOfEdgelx[ln-locVol-bord_vol]);
		const int neighVal=test[ln];
		
		if(neighVal!=gn)
		  master_printf("site %s ori (%d,%d) dir (%d,%d) has edgelx neigh %s with val %s\n",
				siteAsString(glblxOfLoclx[site]).c_str(),
				ori1,ori2,mu,nu,
				siteAsString(gn).c_str(),
				siteAsString(neighVal).c_str());
	      }
      }
    NISSA_PARALLEL_LOOP_END;
    
    master_printf("lx edges communicates consistently\n");
  }
  
  void testEoEdgesExchange()
  {
    EoField<int> test("testEoEdge",WITH_HALO_EDGES);
    forBothParities([&test](const auto& par)
    {
      NISSA_PARALLEL_LOOP(i,0,locVolh)
	test[par][i]=glblxOfLoclx[loclx_of_loceo[par][i]];
      NISSA_PARALLEL_LOOP_END;
      
      NISSA_PARALLEL_LOOP(i,0,bord_volh)
	test[par][i+locVolh]=-1;
      NISSA_PARALLEL_LOOP_END;
      
      NISSA_PARALLEL_LOOP(i,0,edge_volh)
	test[par][i+locVolh+bord_volh]=-2;
      NISSA_PARALLEL_LOOP_END;
    });
    
    int* r=(int*)recv_buf;
    int* s=(int*)send_buf;
    const int n=recv_buf_size/sizeof(int);
    NISSA_PARALLEL_LOOP(i,0,n)
      {
	r[i]=-3;
	s[i]=-4;
      }
    NISSA_PARALLEL_LOOP_END;
    
    test.invalidateHalo();
    test.invalidateEdges();
    test.updateHalo();
    test.updateEdges();
    
    forBothParities([&test](const auto& par)
    {
      NISSA_PARALLEL_LOOP(site,0,locVolh)
	{
	  for(int ori1=0;ori1<2;ori1++)
	    for(int ori2=0;ori2<2;ori2++)
	      for(int iEdge=0;iEdge<nEdges;iEdge++)
		{
		  const auto [mu,nu]=edge_dirs[iEdge];
		  
		  const int l1n=((ori1==0)?loceo_neighdw:loceo_neighup)[par][site][mu];
		  const int ln=((ori2==0)?loceo_neighdw:loceo_neighup)[!par][l1n][nu];
		  const int lx=loclx_of_loceo[par][ln];
		  const int gn=(lx<locVol)?glblxOfLoclx[lx]:((lx<locVol+bord_vol)?glblxOfBordlx[lx-locVol]:glblxOfEdgelx[lx-locVol-bord_vol]);
		  const int neighVal=test[par][ln];
		  
		  if(neighVal!=gn)
		    master_printf("par %d site %s ori (%d,%d) dir (%d,%d) neigh %s with val %s\n",
				  par(),
				  siteAsString(glblxOfLoclx[loclx_of_loceo[par][site]]).c_str(),
				  ori1,ori2,mu,nu,
				  siteAsString(gn).c_str(),
				  siteAsString(neighVal).c_str());
		}
	}
      NISSA_PARALLEL_LOOP_END;
    });
    
    master_printf("eo edges communicates consistently\n");
  }
}
