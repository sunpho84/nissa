#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <string.h>

#include "base/global_variables.hpp"
#include "base/thread_macros.hpp"
#include "base/vectors.hpp"
#include "new_types/complex.hpp"
#include "new_types/dirac.hpp"
#include "new_types/spin.hpp"
#include "new_types/su3.hpp"
#include "routines/ios.hpp"
#ifdef USE_THREADS
 #include "routines/thread.hpp"
#endif

#include "site_contract.hpp"

namespace nissa
{
  //product of the color traces
  THREADABLE_FUNCTION_10ARG(trace_g_sdag_g_s_times_trace_g_sdag_g_s, complex**,glb_c, dirac_matr*,g1L, colorspinspin*,s1L, dirac_matr*,g2L, colorspinspin*,s2L, dirac_matr*,g1R,colorspinspin*,s1R, dirac_matr*,g2R, colorspinspin*,s2R, int,ncontr)
  {
    GET_THREAD_ID();
    
    //allocate a contiguous memory area where to store local results
    complex *loc_c=nissa_malloc("loc_c",ncontr*glb_size[0],complex);
    vector_reset(loc_c);
    
    //loop over time and number of contractions
    NISSA_PARALLEL_LOOP(ibase,0,ncontr*loc_size[0])
      {
	//find out contraction index and time
	int icontr=ibase/loc_size[0];
	int t=ibase-icontr*loc_size[0];
	
	//loop over spatial volume
	for(int ivol=t*loc_spat_vol;ivol<(t+1)*loc_spat_vol;ivol++)
	  {
	    int glb_t=glb_coord_of_loclx[ivol][0];
	    
	    //color loops on the left
	    complex ctempL;
	    trace_g_css_dag_g_css(ctempL,&(g1L[icontr]),s1L[ivol],&(g2L[icontr]),s2L[ivol]);
	    
	    //color loop on the right
	    complex ctempR;
	    trace_g_css_dag_g_css(ctempR,&(g1R[icontr]),s1R[ivol],&(g2R[icontr]),s2R[ivol]);
	    
	    //product of the two
	    complex_summ_the_prod(loc_c[icontr*glb_size[0]+glb_t],ctempL,ctempR);
	  }
      }
    
    //wait that all threads finished their piece
    THREAD_BARRIER();
    
    //final reduction
    if(IS_MASTER_THREAD)
      {
	verbosity_lv3_master_printf("Performing final reduction of %d bytes\n",2*glb_size[0]*ncontr);
	MPI_Reduce(loc_c,glb_c,2*glb_size[0]*ncontr,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	verbosity_lv3_master_printf("Reduction done\n");
      }
    
    nissa_free(loc_c);
  }}

  //trace of the product
  THREADABLE_FUNCTION_10ARG(trace_g_css_dag_g_ss_g_css_dag_g_ss, complex**,glb_c, dirac_matr*,g1L, colorspinspin*,s1L, dirac_matr*,g2L, colorspinspin*,s2L, dirac_matr*,g1R,colorspinspin*,s1R, dirac_matr*,g2R, colorspinspin*,s2R, int,ncontr)
  {
    GET_THREAD_ID();
    
    //allocate a contiguous memory area where to store local results
    complex *loc_c=nissa_malloc("loc_c",ncontr*glb_size[0],complex);
    vector_reset(loc_c);
    
    //loop over time and number of contractions
    NISSA_PARALLEL_LOOP(ibase,0,ncontr*loc_size[0])
      {
	//find out contraction index and time
	int icontr=ibase/loc_size[0];
	int t=ibase-icontr*loc_size[0];
	
	//loop over spatial volume
	for(int ivol=t*loc_spat_vol;ivol<(t+1)*loc_spat_vol;ivol++)
	  {
	    int glb_t=glb_coord_of_loclx[ivol][0];
	    
	    //color loops
	    for(int icol1=0;icol1<3;icol1++)
	      for(int icol2=0;icol2<3;icol2++)
		{
		  complex ctemp;
		  trace_g_ss_dag_g_ss_g_ss_dag_g_ss(ctemp,&(g1L[icontr]),s1L[ivol][icol2],&(g2L[icontr]),s2R[ivol][icol2],
						    &(g1R[icontr]),s1R[ivol][icol1],&(g2R[icontr]),s2L[ivol][icol1]);
		  complex_summassign(loc_c[icontr*glb_size[0]+glb_t],ctemp);
		}
	  }
      }
    
    //wait that all threads finished their part
    THREAD_BARRIER();
    
    //final reduction
    if(IS_MASTER_THREAD)
      {
	verbosity_lv3_master_printf("Performing final reduction of %d bytes\n",2*glb_size[0]*ncontr);
	MPI_Reduce(loc_c,glb_c,2*glb_size[0]*ncontr,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	verbosity_lv3_master_printf("Reduction done\n");
      }
    
    nissa_free(loc_c);
  }}

  THREADABLE_FUNCTION_8ARG(trace_id_css_dag_g_css_id_css_dag_g_css, complex*,glb_c, colorspinspin*,s1L, dirac_matr*,g2L, colorspinspin*,s2L, colorspinspin*,s1R, dirac_matr*,g2R, colorspinspin*,s2R, int,ncontr)
  {
    GET_THREAD_ID();
    
    //allocate a contiguous memory area where to store local results
    complex *loc_c=nissa_malloc("loc_c",ncontr*glb_size[0],complex);
    vector_reset(loc_c);
    
    //loop over time and number of contractions
    NISSA_PARALLEL_LOOP(ibase,0,ncontr*loc_size[0])
      {
	//find out contraction index and time
	int icontr=ibase/loc_size[0];
	int t=ibase-icontr*loc_size[0];
	int glb_t=glb_coord_of_loclx[0][0]+t;
	
	//loop over spatial volume
	for(int ivol=t*loc_spat_vol;ivol<(t+1)*loc_spat_vol;ivol++)
	  for(int icol1=0;icol1<3;icol1++)
	    for(int icol2=0;icol2<3;icol2++)
	      {
		spinspin AR,AL;
		unsafe_spinspin_prod_spinspin_dag(AL,s2L[ivol][icol1],s1L[ivol][icol2]);
		unsafe_spinspin_prod_spinspin_dag(AR,s2R[ivol][icol2],s1R[ivol][icol1]);
		
		for(int icontr=0;icontr<ncontr;icontr++)
		  {
		    complex ctemp;
		    spinspin ALg, ARg;
		    
		    unsafe_dirac_prod_spinspin(ALg,g2R+icontr,AL);
		    unsafe_dirac_prod_spinspin(ARg,g2L+icontr,AR);
		    trace_prod_spinspins(ctemp,ALg,ARg);
		    complex_summassign(loc_c[icontr*glb_size[0]+glb_t],ctemp);
		  }
	      }
      }
    
    THREAD_BARRIER();
    
    //final reduction
    if(IS_MASTER_THREAD)
      {
	verbosity_lv3_master_printf("Performing final reduction of %d bytes\n",2*glb_size[0]*ncontr);
	MPI_Reduce(loc_c,glb_c,2*glb_size[0]*ncontr,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	verbosity_lv3_master_printf("Reduction done\n");
      }
    
    nissa_free(loc_c);
  }}

  THREADABLE_FUNCTION_8ARG(trace_id_css_dag_g_css_times_trace_id_css_dag_g_css, complex*,glb_c, colorspinspin*,s1L, dirac_matr*,g2L, colorspinspin*,s2L, colorspinspin*,s1R, dirac_matr*,g2R, colorspinspin*,s2R, int,ncontr)
  {
    GET_THREAD_ID();
    
    //allocate a contiguous memory area where to store local results
    complex *loc_c=nissa_malloc("loc_c",ncontr*glb_size[0],complex);
    vector_reset(loc_c);
    
    //loop over time and number of contractions
    NISSA_PARALLEL_LOOP(ibase,0,ncontr*loc_size[0])
      {
	//find out contraction index and time
	int icontr=ibase/loc_size[0];
	int t=ibase-icontr*loc_size[0];
	
	//loop over spatial volume
	for(int ivol=t*loc_spat_vol;ivol<(t+1)*loc_spat_vol;ivol++)
	  {
	    int glb_t=glb_coord_of_loclx[ivol][0];
	    complex ctempL[ncontr];
	    complex ctempR[ncontr];
	    complex ctemp[ncontr];
	    memset(ctempL,0,ncontr*sizeof(complex));
	    memset(ctempR,0,ncontr*sizeof(complex));
	    memset(ctemp,0,ncontr*sizeof(complex));
	    
	    for(int icol=0;icol<3;icol++)
	      {
		spinspin AL,AR;
		unsafe_spinspin_prod_spinspin_dag(AL,s2L[ivol][icol],s1L[ivol][icol]);
		unsafe_spinspin_prod_spinspin_dag(AR,s2R[ivol][icol],s1R[ivol][icol]);
		for(int icontr=0;icontr<ncontr;icontr++)
		  {
		    complex ctempL_color,ctempR_color;
		    
		    trace_dirac_prod_spinspin(ctempL_color,g2R+icontr,AL);
		    trace_dirac_prod_spinspin(ctempR_color,g2L+icontr,AR);
		    
		    complex_summ(ctempL[icontr],ctempL[icontr],ctempL_color);
		    complex_summ(ctempR[icontr],ctempR[icontr],ctempR_color);
		  }
	      }
	    
	    for(int icontr=0;icontr<ncontr;icontr++)
	      complex_summ_the_prod(loc_c[icontr*glb_size[0]+glb_t],ctempL[icontr],ctempR[icontr]);
	  }
      }
    
    THREAD_BARRIER();
    
    //Final reduction
    if(IS_MASTER_THREAD)
      {
	verbosity_lv3_master_printf("Performing final reduction of %d bytes\n",2*glb_size[0]*ncontr);
	MPI_Reduce(loc_c,glb_c,2*glb_size[0]*ncontr,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	verbosity_lv3_master_printf("Reduction done\n");
      }
    
    nissa_free(loc_c);
  }}
}
