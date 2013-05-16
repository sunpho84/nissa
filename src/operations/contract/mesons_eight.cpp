#ifdef HAVE_CONFIG_H
 #include "config.h"
#endif

#include <math.h>
#include <string.h>

#include "../base/global_variables.h"
#include "../base/vectors.h"
#include "../communicate/communicate.h"
#include "../new_types/complex.h"
#include "../new_types/dirac.h"
#include "../new_types/new_types_definitions.h"
#include "../new_types/spin.h"
#include "../new_types/su3.h"
#include "../routines/ios.h"
#include "../routines/thread.h"

//the eight
void site_trace_g_sdag_g_s_g_sdag_g_s(complex c,dirac_matr *g1,spinspin s1,dirac_matr *g2,spinspin s2, dirac_matr *g3,spinspin s3,dirac_matr *g4,spinspin s4)
{
  spinspin t1,t2,t12,t3,t4,t34;
  
  unsafe_dirac_prod_spinspin_dag(t1,g1,s1);
  unsafe_dirac_prod_spinspin(t2,g2,s2);
  unsafe_spinspin_prod_spinspin(t12,t1,t2);

  unsafe_dirac_prod_spinspin_dag(t3,g3,s3);
  unsafe_dirac_prod_spinspin(t4,g4,s4);
  unsafe_spinspin_prod_spinspin(t34,t3,t4);
  
  trace_prod_spinspins(c,t12,t34);
}

void sum_trace_g_sdag_g_s_times_trace_g_sdag_g_s(complex **glb_c, dirac_matr *g1L,colorspinspin *s1L, dirac_matr *g2L, colorspinspin *s2L, dirac_matr *g1R,colorspinspin *s1R, dirac_matr *g2R, colorspinspin *s2R,const int ncontr)
{
//Allocate a contiguous memory area where to store local results
  complex *loc_c=nissa_malloc("loc_c",ncontr*glb_size[0],complex);
  for(int icontr=0;icontr<ncontr;icontr++)
    for(int glb_t=0;glb_t<glb_size[0];glb_t++) loc_c[icontr*glb_size[0]+glb_t][0]=loc_c[icontr*glb_size[0]+glb_t][1]=0;

  //Check if the global vector is contiguos
  complex *glb_c_buf=glb_c[0];
  int use_buf=0;
  for(int icontr=0;icontr<ncontr && use_buf==0;icontr++) use_buf=(glb_c[icontr]!=glb_c_buf+icontr*glb_size[0]);
  if(use_buf)
    {
      verbosity_lv2_master_printf("Creating a temporary buffer for the contractions.\n");
      verbosity_lv2_master_printf("Avoid it passing a 'glb_c' pointing to a contiguos memory area\n");
      glb_c_buf=nissa_malloc("glb_c_buf",ncontr*glb_size[0],complex);
    }

  for(int icontr=0;icontr<ncontr;icontr++)
    {
      verbosity_lv3_master_printf("Contraction %d/%d\n",icontr+1,ncontr);

      //Local loop
      nissa_loc_vol_loop(ivol)
        {
          complex ctemp;
	  complex ctempL_color,ctempL={0,0};
	  complex ctempR_color,ctempR={0,0};
	  
          int glb_t=glb_coord_of_loclx[ivol][0];
          //Color loops
          for(int icol=0;icol<3;icol++)
	    {
	      site_trace_g_sdag_g_s(ctempL_color,&(g1L[icontr]),s1L[ivol][icol],&(g2L[icontr]),s2L[ivol][icol]);
	      complex_summ(ctempL,ctempL,ctempL_color);
	    }
	 for(int icol=0;icol<3;icol++)
	   {
	     site_trace_g_sdag_g_s(ctempR_color,&(g1R[icontr]),s1R[ivol][icol],&(g2R[icontr]),s2R[ivol][icol]);
	     complex_summ(ctempR,ctempR,ctempR_color);
	   }
	 safe_complex_prod(ctemp,ctempL,ctempR);
	 complex_summ(loc_c[icontr*glb_size[0]+glb_t],loc_c[icontr*glb_size[0]+glb_t],ctemp);
        }
    }

  //final reduction
  verbosity_lv3_master_printf("Performing final reduction of %d bytes\n",2*glb_size[0]*ncontr);
  MPI_Reduce(loc_c,glb_c_buf,2*glb_size[0]*ncontr,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
  verbosity_lv3_master_printf("Reduction done\n");
  
  //if a temporary buffer has been used, destory it after copyng data to the true one
  if(use_buf)
    {
      for(int icontr=0;icontr<ncontr;icontr++)
        for(int glb_t=0;glb_t<glb_size[0];glb_t++)
          {
            glb_c[icontr][glb_t][0]=glb_c_buf[icontr*glb_size[0]+glb_t][0];
            glb_c[icontr][glb_t][1]=glb_c_buf[icontr*glb_size[0]+glb_t][1];
          }
      nissa_free(glb_c_buf);
    }

  nissa_free(loc_c);
}

void trace_g_sdag_g_s_g_sdag_g_s(complex **glb_c, dirac_matr *g1L,colorspinspin *s1L, dirac_matr *g2L, colorspinspin *s2L, dirac_matr *g1R,colorspinspin *s1R, dirac_matr *g2R, colorspinspin *s2R,const int ncontr)
{
//Allocate a contguous memory area where to store local results
  complex *loc_c=nissa_malloc("loc_c",ncontr*glb_size[0],complex);
  for(int icontr=0;icontr<ncontr;icontr++)
    for(int glb_t=0;glb_t<glb_size[0];glb_t++) loc_c[icontr*glb_size[0]+glb_t][0]=loc_c[icontr*glb_size[0]+glb_t][1]=0;

  //Check if the global vector is contiguos
  complex *glb_c_buf=glb_c[0];
  int use_buf=0;
  for(int icontr=0;icontr<ncontr && use_buf==0;icontr++) use_buf=(glb_c[icontr]!=glb_c_buf+icontr*glb_size[0]);
  if(use_buf)
    {
      verbosity_lv1_master_printf("Creating a temporary buffer for the contractions.\n");
      verbosity_lv1_master_printf("Avoid it passing a 'glb_c' pointing to a contiguos memory area\n");
      glb_c_buf=nissa_malloc("glb_c_buf",ncontr*glb_size[0],complex);
    }

  for(int icontr=0;icontr<ncontr;icontr++)
    {
      verbosity_lv3_master_printf("Contraction %d/%d\n",icontr+1,ncontr);

      //Local loop
      nissa_loc_vol_loop(ivol)
        {
          complex ctemp;
          int glb_t=glb_coord_of_loclx[ivol][0];
          //Color loop
          for(int icol1=0;icol1<3;icol1++) for(int icol2=0;icol2<3;icol2++)
            {
//              site_trace_g_sdag_g_s_g_sdag_g_s(ctemp,&(g1L[icontr]),s1L[ivol][icol1],&(g2L[icontr]),s2L[ivol][icol2],&(g1R[icontr]),s1R[ivol][icol2],&(g2R[icontr]),s2R[ivol][icol1]);
		site_trace_g_sdag_g_s_g_sdag_g_s(ctemp,&(g1L[icontr]),s1L[ivol][icol2],&(g2L[icontr]),s2R[ivol][icol2],&(g1R[icontr]),s1R[ivol][icol1],&(g2R[icontr]),s2L[ivol][icol1]);

              complex_summ(loc_c[icontr*glb_size[0]+glb_t],loc_c[icontr*glb_size[0]+glb_t],ctemp);
            }
        }
    }

  //Final reduction
  verbosity_lv3_master_printf("Performing final reduction of %d bytes\n",2*glb_size[0]*ncontr);
  MPI_Reduce(loc_c,glb_c_buf,2*glb_size[0]*ncontr,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
  verbosity_lv3_master_printf("Reduction done\n");

  //if a temporary buffer has been used, destory it after copyng data to the true one
  if(use_buf)
    {
      for(int icontr=0;icontr<ncontr;icontr++)
        for(int glb_t=0;glb_t<glb_size[0];glb_t++)
          {
            glb_c[icontr][glb_t][0]=glb_c_buf[icontr*glb_size[0]+glb_t][0];
            glb_c[icontr][glb_t][1]=glb_c_buf[icontr*glb_size[0]+glb_t][1];
          }
      nissa_free(glb_c_buf);
    }

  nissa_free(loc_c);
}

void trace_id_sdag_g_s_id_sdag_g_s(complex *glb_c,colorspinspin *s1L,dirac_matr *g2L,colorspinspin *s2L,colorspinspin *s1R,dirac_matr *g2R,colorspinspin *s2R,const int ncontr)
{
  //Allocate a contiguous memory area where to store local results
  complex *loc_c=nissa_malloc("loc_c",ncontr*glb_size[0],complex);
  vector_reset(loc_c);

  //Local loop
  spinspin AR,AL;
  nissa_loc_vol_loop(ivol)
    {
      int glb_t=glb_coord_of_loclx[ivol][0];
      for(int icol1=0;icol1<3;icol1++)
	for(int icol2=0;icol2<3;icol2++)
	  {
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

  //Final reduction
  verbosity_lv3_master_printf("Performing final reduction of %d bytes\n",2*glb_size[0]*ncontr);
  MPI_Reduce(loc_c,glb_c,2*glb_size[0]*ncontr,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
  verbosity_lv3_master_printf("Reduction done\n");
  
  nissa_free(loc_c);
}

void sum_trace_id_sdag_g_s_times_trace_id_sdag_g_s(complex *glb_c,colorspinspin *s1L,dirac_matr *g2L,colorspinspin *s2L,colorspinspin *s1R,dirac_matr *g2R,colorspinspin *s2R,const int ncontr)
{
  //Allocate a contguous memory area where to store local results
  complex *loc_c=nissa_malloc("loc_c",ncontr*glb_size[0],complex);
  vector_reset(loc_c);
  
  //Local loop
  nissa_loc_vol_loop(ivol)
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
	      trace_dirac_prod_spinspin(ctempR_color,&(g2L[icontr]),AR);

	      complex_summ(ctempL[icontr],ctempL[icontr],ctempL_color);
	      complex_summ(ctempR[icontr],ctempR[icontr],ctempR_color);
	    }
	}
      
      for(int icontr=0;icontr<ncontr;icontr++)
	complex_summ_the_prod(loc_c[icontr*glb_size[0]+glb_t],ctempL[icontr],ctempR[icontr]);
    }
  
  //Final reduction
  verbosity_lv3_master_printf("Performing final reduction of %d bytes\n",2*glb_size[0]*ncontr);
  MPI_Reduce(loc_c,glb_c,2*glb_size[0]*ncontr,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
  verbosity_lv3_master_printf("Reduction done\n");
  
  nissa_free(loc_c);
}
