#ifdef HAVE_CONFIG_H
 #include "config.h"
#endif

#include <math.h>
#include <string.h>

#include "../../base/global_variables.h"
#include "../../base/vectors.h"
#include "../../new_types/complex.h"
#include "../../new_types/dirac.h"
#include "../../new_types/su3.h"
#include "../../routines/ios.h"
#include "../../routines/thread.h"

//take Tr[g1 * s1^dag * g2 * s2], useful for mesons 2 points
void trace_g_ss_dag_g_ss(complex c,dirac_matr *g1,spinspin s1,dirac_matr *g2,spinspin s2)
{
  spinspin t1,t2;
  
  unsafe_dirac_prod_spinspin_dag(t1,g1,s1);
  unsafe_dirac_prod_spinspin(t2,g2,s2);

  trace_prod_spinspins(c,t1,t2);
}
void trace_g_css_dag_g_css(complex c,dirac_matr *g1,colorspinspin s1,dirac_matr *g2,colorspinspin s2)
{
  //reset out
  c[0]=c[1]=0;
  
  //loop over color indices
  for(int ic=0;ic<3;ic++)
    {
      spinspin t1,t2;
      
      unsafe_dirac_prod_spinspin_dag(t1,g1,s1[ic]);
      unsafe_dirac_prod_spinspin(t2,g2,s2[ic]);
      
      summ_the_trace_prod_spinspins(c,t1,t2);
    }
}
void trace_g_ccss_dag_g_ccss(complex c,dirac_matr *g1,su3spinspin s1,dirac_matr *g2,su3spinspin s2)
{
  //reset out
  c[0]=c[1]=0;
  
  //loop over color indices
  for(int ic1=0;ic1<3;ic1++)
    for(int ic2=0;ic2<3;ic2++)
      {
	spinspin t1,t2;
	
	unsafe_dirac_prod_spinspin_dag(t1,g1,s1[ic2][ic1]);
	unsafe_dirac_prod_spinspin(t2,g2,s2[ic2][ic1]);
	
	summ_the_trace_prod_spinspins(c,t1,t2);
      }
}

#define DEFINE_TWO_POINTS_MESON_ROUTINES_FOR_TYPE(TYPE,SHORTTYPE)	\
  THREADABLE_FUNCTION_6ARG(NAME4(trace_g,SHORTTYPE,dag_g,SHORTTYPE), complex*,glb_c, dirac_matr*,g1, TYPE*,s1, dirac_matr*,g2, TYPE*,s2, int,ncontr) \
  {									\
    GET_THREAD_ID();							\
    									\
    /*local buffer*/							\
    complex *loc_c=nissa_malloc("loc_c",glb_size[0]*ncontr,complex);	\
									\
    /*loop over time and number of contractions*/			\
    NISSA_PARALLEL_LOOP(ibase,0,ncontr*loc_size[0])			\
      {									\
	/*find out contraction index and time*/				\
	int icontr=ibase/loc_size[0];					\
	int t=ibase-icontr*loc_size[0];					\
									\
	/*loop over spatial volume*/					\
	for(int ivol=t*loc_spat_vol;ivol<(t+1)*loc_spat_vol;ivol++)	\
	  {								\
	    int glb_t=glb_coord_of_loclx[ivol][0];			\
									\
	    /*summ the trace*/						\
	    complex ctemp;						\
	    NAME4(trace_g,SHORTTYPE,dag_g,SHORTTYPE)(ctemp,&(g1[icontr]),s1[ivol],&(g2[icontr]),s2[ivol]); \
	    complex_summassign(loc_c[icontr*glb_size[0]+glb_t],ctemp);	\
	  }								\
      }									\
									\
    /*wait that all threads finish*/					\
    THREAD_BARRIER();							\
									\
    if(IS_MASTER_THREAD)						\
      {									\
	verbosity_lv3_master_printf("Performing final reduction of %d double\n",2*glb_size[0]*ncontr); \
	MPI_Reduce(loc_c,glb_c,2*glb_size[0]*ncontr,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD); \
	verbosity_lv3_master_printf("Reduction done\n");		\
      }									\
									\
    /*free local buffer and auto sync*/					\
   nissa_free(loc_c);							\
  }}

//this function takes care to make the revert on the FIRST spinor, putting the needed gamma5
#define DEFINE_MESON_TWO_POINTS_WILSON_PROP(TYPE,SHORTTYPE)		\
  void meson_two_points_Wilson_prop(complex *corr,int *list_op_source,TYPE *s1,int *list_op_sink,TYPE *s2,int ncontr) \
  {									\
    /*temporary vectors for the internal gamma*/			\
    dirac_matr tsource[ncontr],tsink[ncontr];				\
    									\
    for(int icontr=0;icontr<ncontr;icontr++)				\
      {									\
	/*put the two gamma5 needed for the revert of the first spinor*/ \
	dirac_prod(&(tsource[icontr]), &(base_gamma[list_op_source[icontr]]),&(base_gamma[5]));	\
	dirac_prod(&(tsink[icontr]), &(base_gamma[5]),&(base_gamma[list_op_sink[icontr]])); \
      }									\
									\
    /*call the routine which perform the contraction*/			\
    NAME4(trace_g,SHORTTYPE,dag_g,SHORTTYPE)(corr,tsource,s1,tsink,s2,ncontr); \
}

DEFINE_TWO_POINTS_MESON_ROUTINES_FOR_TYPE(colorspinspin,css)
DEFINE_TWO_POINTS_MESON_ROUTINES_FOR_TYPE(su3spinspin,ccss)

DEFINE_MESON_TWO_POINTS_WILSON_PROP(colorspinspin,css)
DEFINE_MESON_TWO_POINTS_WILSON_PROP(su3spinspin,ccss)
