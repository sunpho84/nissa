#ifdef HAVE_CONFIG_H
# include "config.hpp"
#endif

#include <math.h>
#include <string.h>
#include <fftw3.h>

#include "base/debug.hpp"
#include "base/vectors.hpp"
#include "geometry/geometry_lx.hpp"
#include "linalgs/linalgs.hpp"
#include "operations/remap_vector.hpp"
#include "routines/ios.hpp"
#include "routines/math_routines.hpp"
#include "routines/mpi_routines.hpp"

namespace nissa
{
  void fft4d(complex *out,
	     const complex *in,
	     const WhichDirs& dirs,
	     const int& ncpp,
	     const double& sign,
	     const bool& normalize)
  {
    //first of all put in to out
    if(out!=in) vector_copy(out,in);
    
    //list all dirs
    int list_dirs[NDIM],ndirs=0;
    for(int mu=0;mu<NDIM;mu++) if(dirs[mu]) list_dirs[ndirs++]=mu;
    VERBOSITY_LV2_MASTER_PRINTF("Going to FFT: %d dimensions in total\n",ndirs);
    
    if(ndirs)
      {
	//allocate buffer
	complex *buf=nissa_malloc("buf",max_locd_size*ncpp,complex);
	
	//allocate plans
	fftw_plan *plans=nissa_malloc("plans",ndirs,fftw_plan);
	  for(int idir=0;idir<ndirs;idir++)
	    plans[idir]=fftw_plan_many_dft(1,&glbSize[list_dirs[idir]],ncpp,buf,NULL,ncpp,1,buf,NULL,ncpp,1,sign,FFTW_ESTIMATE);
	
	//transpose each dir in turn and take fft
	for(int idir=0;idir<ndirs;idir++)
	  {
	    int mu=list_dirs[idir];
	    VERBOSITY_LV2_MASTER_PRINTF("FFT-ing dimension %d/%d=%d\n",idir+1,ndirs,mu);
	    remap_lx_vector_to_locd(buf,out,ncpp*sizeof(complex),mu);
	    
	    //makes all the fourier transform
	    HOST_PARALLEL_LOOP(0,locd_perp_size_per_dir[mu],
			       CAPTURE(plans,idir,buf,mu,ncpp),ioff,
				{
				  fftw_execute_dft(plans[idir],buf+ioff*glbSize[mu]*ncpp,buf+ioff*glbSize[mu]*ncpp);
				});
	    
	    remap_locd_vector_to_lx(out,buf,ncpp*sizeof(complex),mu);
	  }
	
	//destroy plans
	for(int idir=0;idir<ndirs;idir++) fftw_destroy_plan(plans[idir]);
	
	//put normaliisation
	if(normalize)
	  {
	    double norm=glbSize[list_dirs[0]];
	    for(int idir=1;idir<ndirs;idir++) norm*=glbSize[idir];
	    PAR(0,locVol*ncpp,
		CAPTURE(out,norm),
		i,
		{
		  complex_prodassign_double(out[i],1/norm);
		});
	  }
	
	nissa_free(buf);
	nissa_free(plans);
      }
  }
}
