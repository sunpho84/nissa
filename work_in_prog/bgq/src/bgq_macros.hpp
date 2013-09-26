#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "new_vars_and_types.hpp"

////////////////////////////////////////////// debug ////////////////////////////////////////////

#include <math.hpp>

#define CHECK_SU3_DIFF_FROM_ZERO(U)					
  /*
  {                                                                     \
    for(int ic1=0;ic1<3;ic1++)						\
      {									\
	complex d={0,0};						\
	for(int ic2=0;ic2<3;ic2++)					\
	  complex_summ_thppe_conj2_prod(d,U[ic1][ic2],U[ic1][ic2]); \
	if(fabs(d[0]-1)>1.e-14)						\
	  crashpp("su3 not unitary, vn %d, ic1 %d, d=%lg!",ic1,d[0]);	\
      }									\
  }
  */
#define CHECK_BI_SU3_DIFF_FROM_ZERO(U)					
/*
  {									\
    for(int vn=0;vn<2;vn++)						\
      for(int ic1=0;ic1<3;ic1++)					\
	{								\
	  complex d={0,0};						\
	  for(int ic2=0;ic2<3;ic2++)					\
	    complex_summ_thppe_conj2_prod(d,U[ic1][ic2][vn],U[ic1][ic2][vn]); \
	  if(fabs(d[0]-1)>1.e-14)					\
	    crashpp("not unitary, vn %d, ic1 %d, d=%lg!",vn,ic1,d[0]);	\
	}								\
  }
*/

#define CHECK_BI_SPINCOLOR_DIFF_FROM_ZERO(S,ivol)			
/*
  {									\
    for(int id=0;id<2;id++)						\
      for(int ic=0;ic<3;ic++)						\
	for(int vn=0;vn<2;vn++)						\
	  for(int ri=0;ri<2;ri++)					\
	    if(fabs(S[id][ic][vn][ri])<1.e-14)				\
	      crashpp("ivol: %d id %d ic: %d vn: %d ri: %d is empty (%lg)",ivol,id,ic,vn,ri,out[ivol][id][ic][vn][ri]); \
  }
*/
////////////////////////////// convert normal complex to bi_complex //////////////////////////

#define COMPLEX_TO_BI_COMPLEX(A,B,VN)		\
  {						\
    complex_copy(A[VN],B);			\
  }
#define COLOR_TO_BI_COLOR(A,B,VN)		\
  {						\
    COMPLEX_TO_BI_COMPLEX(A[0],B[0],VN);	\
    COMPLEX_TO_BI_COMPLEX(A[1],B[1],VN);	\
    COMPLEX_TO_BI_COMPLEX(A[2],B[2],VN);	\
  }
#define SU3_TO_BI_SU3(A,B,VN)			\
  {						\
    CHECK_SU3_DIFF_FROM_ZERO(B);		\
    COLOR_TO_BI_COLOR(A[0],B[0],VN);		\
    COLOR_TO_BI_COLOR(A[1],B[1],VN);		\
    COLOR_TO_BI_COLOR(A[2],B[2],VN);		\
  }
#define HALFSPINCOLOR_TO_BI_HALFSPINCOLOR(A,B,VN)	\
  {							\
    COLOR_TO_BI_COLOR(A[0],B[0],VN);			\
    COLOR_TO_BI_COLOR(A[1],B[1],VN);			\
  }
#define SPINCOLOR_TO_BI_SPINCOLOR(A,B,VN)	\
  {						\
    COLOR_TO_BI_COLOR(A[0],B[0],VN);		\
    COLOR_TO_BI_COLOR(A[1],B[1],VN);		\
    COLOR_TO_BI_COLOR(A[2],B[2],VN);		\
    COLOR_TO_BI_COLOR(A[3],B[3],VN);		\
  }

//////////////////////////////// extract component of bi_complex ////////////////////////////////

#define COMPLEX_OF_BI_COMPLEX(A,B,VN)		\
  {						\
    complex_copy(A,B[VN]);			\
  }
#define COLOR_OF_BI_COLOR(A,B,VN)		\
  {						\
    COMPLEX_OF_BI_COMPLEX(A[0],B[0],VN);	\
    COMPLEX_OF_BI_COMPLEX(A[1],B[1],VN);	\
    COMPLEX_OF_BI_COMPLEX(A[2],B[2],VN);	\
  }
#define SU3_OF_BI_SU3(A,B,VN)			\
  {						\
    CHECK_SU3_DIFF_FROM_ZERO(B);		\
    COLOR_OF_BI_COLOR(A[0],B[0],VN);		\
    COLOR_OF_BI_COLOR(A[1],B[1],VN);		\
    COLOR_OF_BI_COLOR(A[2],B[2],VN);		\
  }
#define HALFSPINCOLOR_OF_BI_HALFSPINCOLOR(A,B,VN)	\
  {							\
    COLOR_OF_BI_COLOR(A[0],B[0],VN);			\
    COLOR_OF_BI_COLOR(A[1],B[1],VN);			\
  }
#define SPINCOLOR_OF_BI_SPINCOLOR(A,B,VN)	\
  {						\
    COLOR_OF_BI_COLOR(A[0],B[0],VN);		\
    COLOR_OF_BI_COLOR(A[1],B[1],VN);		\
    COLOR_OF_BI_COLOR(A[2],B[2],VN);		\
    COLOR_OF_BI_COLOR(A[3],B[3],VN);		\
  }

//////////////////////////////// copy bi_complex ////////////////////////////////

#define BI_COMPLEX_COPY(A,B)			\
  {						\
    complex_copy(A[0],B[0]);			\
    complex_copy(A[1],B[1]);			\
  }
#define BI_COLOR_COPY(A,B)			\
  {						\
    BI_COMPLEX_COPY(A[0],B[0]);			\
    BI_COMPLEX_COPY(A[1],B[1]);			\
    BI_COMPLEX_COPY(A[2],B[2]);			\
  }
#define BI_SU3_COPY(A,B)			\
  {						\
    BI_COLOR_COPY(A[0],B[0]);			\
    BI_COLOR_COPY(A[1],B[1]);			\
    BI_COLOR_COPY(A[2],B[2]);			\
  }
#define BI_HALFSPINCOLOR_COPY(A,B)		\
  {						\
    BI_COLOR_COPY(A[0],B[0]);			\
    BI_COLOR_COPY(A[1],B[1]);			\
  }
#define BI_SPINCOLOR_COPY(A,B)			\
  {						\
    BI_COLOR_COPY(A[0],B[0]);			\
    BI_COLOR_COPY(A[1],B[1]);			\
    BI_COLOR_COPY(A[2],B[2]);			\
    BI_COLOR_COPY(A[3],B[3]);			\
  }

/////////////////////////////////// bi_complex to normal type ///////////////////////////////////

#define BI_COMPLEX_TO_COMPLEX(A,B,C)		\
  {						\
    complex_copy(A,C[0]);			\
    complex_copy(B,C[1]);			\
  }
#define BI_COLOR_TO_COLOR(A,B,C)		\
  {						\
    BI_COMPLEX_TO_COMPLEX(A[0],B[0],C[0]);	\
    BI_COMPLEX_TO_COMPLEX(A[1],B[1],C[1]);	\
    BI_COMPLEX_TO_COMPLEX(A[2],B[2],C[2]);	\
  }
#define BI_SU3_TO_SU3(A,B,C)			\
  {						\
    BI_COLOR_TO_COLOR(A[0],B[0],C[0]);		\
    BI_COLOR_TO_COLOR(A[1],B[1],C[1]);		\
    BI_COLOR_TO_COLOR(A[2],B[2],C[2]);		\
  }
#define BI_HALFSPINCOLOR_TO_HALFSPINCOLOR(A,B,C)	\
  {							\
    BI_COLOR_TO_COLOR(A[0],B[0],C[0]);			\
    BI_COLOR_TO_COLOR(A[1],B[1],C[1]);			\
  }
#define BI_SPINCOLOR_TO_SPINCOLOR(A,B,C)	\
  {						\
    BI_COLOR_TO_COLOR(A[0],B[0],C[0]);		\
    BI_COLOR_TO_COLOR(A[1],B[1],C[1]);		\
    BI_COLOR_TO_COLOR(A[2],B[2],C[2]);		\
    BI_COLOR_TO_COLOR(A[3],B[3],C[3]);		\
  }

//////////////////////////////////// basic operation on BI types ////////////////////////////////

#define BI_COMPLEX_PRINT(A)			\
  {						\
    complex_print(A[0]);			\
    complex_print(A[1]);			\
  }

#define BI_COMPLEX_SUBT(A,B,C)                    \
  {						  \
    complex_subt(A[0],B[0],C[0]);		  \
    complex_subt(A[1],B[1],C[1]);		  \
  }
#define BI_COMPLEX_SUMM(A,B,C)                    \
  {						  \
    complex_summ(A[0],B[0],C[0]);		  \
    complex_summ(A[1],B[1],C[1]);		  \
  }

#define BI_COMPLEX_ISUBT(A,B,C)			  \
  {						  \
    complex_isubt(A[0],B[0],C[0]);		  \
    complex_isubt(A[1],B[1],C[1]);		  \
  }

#define BI_COMPLEX_ISUMM(A,B,C)			  \
  {						  \
    complex_isumm(A[0],B[0],C[0]);		  \
    complex_isumm(A[1],B[1],C[1]);		  \
  }

#define BI_COLOR_SUBT(A,B,C)			\
  {						\
    BI_COMPLEX_SUBT(A[0],B[0],C[0]);		\
    BI_COMPLEX_SUBT(A[1],B[1],C[1]);		\
    BI_COMPLEX_SUBT(A[2],B[2],C[2]);		\
  }

#define BI_COLOR_SUMM(A,B,C)			\
  {						\
    BI_COMPLEX_SUMM(A[0],B[0],C[0]);		\
    BI_COMPLEX_SUMM(A[1],B[1],C[1]);		\
    BI_COMPLEX_SUMM(A[2],B[2],C[2]);		\
  }

#define BI_COLOR_ISUBT(A,B,C)			\
  {						\
    BI_COMPLEX_ISUBT(A[0],B[0],C[0]);		\
    BI_COMPLEX_ISUBT(A[1],B[1],C[1]);		\
    BI_COMPLEX_ISUBT(A[2],B[2],C[2]);		\
  }

#define BI_COLOR_ISUMM(A,B,C)			\
  {						\
    BI_COMPLEX_ISUMM(A[0],B[0],C[0]);		\
    BI_COMPLEX_ISUMM(A[1],B[1],C[1]);		\
    BI_COMPLEX_ISUMM(A[2],B[2],C[2]);		\
  }

#define BI_COLOR_PROD_COMPLEX(A,B,C)		\
  {						\
    BI_COMPLEX_PROD(A[0],B[0],C);		\
    BI_COMPLEX_PROD(A[1],B[1],C);		\
    BI_COMPLEX_PROD(A[2],B[2],C);		\
  }

#define BI_COLOR_PROD_DOUBLE(A,B,C)		\
  {						\
    BI_COMPLEX_PROD_DOUBLE(A[0],B[0],C);	\
    BI_COMPLEX_PROD_DOUBLE(A[1],B[1],C);	\
    BI_COMPLEX_PROD_DOUBLE(A[2],B[2],C);	\
  }

#define BI_COLOR_SUMMASSIGN(A,B) BI_COLOR_SUMM(A,A,B)
#define BI_COLOR_SUBTASSIGN(A,B) BI_COLOR_SUBT(A,A,B)
#define BI_COLOR_ISUMMASSIGN(A,B) BI_COLOR_ISUMM(A,A,B)
#define BI_COLOR_ISUBTASSIGN(A,B) BI_COLOR_ISUBT(A,A,B)

//////////////////////////////////////////////////////////////////////////////////////////////////////

#define BI_HALFSPINCOLOR_SUMM(A,B,C)		\
  {						\
    BI_COLOR_SUMM(A[0],B[0],C[0]);		\
    BI_COLOR_SUMM(A[1],B[1],C[1]);		\
  }

#define BI_HALFSPINCOLOR_SUMMASSIGN(A,B)	\
  BI_HALFSPINCOLOR_SUMM(A,A,B);

#define BI_SPINCOLOR_PROD_DOUBLE(A,B,C)					\
  {									\
    BI_COLOR_PROD_DOUBLE(A[0],B[0],C);					\
    BI_COLOR_PROD_DOUBLE(A[1],B[1],C);					\
    BI_COLOR_PROD_DOUBLE(A[2],B[2],C);					\
    BI_COLOR_PROD_DOUBLE(A[3],B[3],C);					\
  }

//////////////////////////////////////////////////////////////////////////////////////////////////////

#define BI_COMPLEX_CONJ1_PROD(A,B,C);		\
  {						\
    unsafe_complex_conj1_prod(A[0],B[0],C[0]);	\
    unsafe_complex_conj1_prod(A[1],B[1],C[1]);	\
  }

#define BI_COMPLEX_PROD(A,B,C);			\
  {						\
    unsafe_complex_prod(A[0],B[0],C[0]);	\
    unsafe_complex_prod(A[1],B[1],C[1]);	\
  }

#define BI_COMPLEX_PROD_DOUBLE(A,B,C);		\
  {						\
    complex_prod_double(A[0],B[0],C);		\
    complex_prod_double(A[1],B[1],C);		\
  }

#define BI_COMPLEX_SUMM_THE_CONJ1_PROD(A,B,C);		\
  {							\
    complex_summ_thppe_conj1_prod(A[0],B[0],C[0]);	\
    complex_summ_thppe_conj1_prod(A[1],B[1],C[1]);	\
  }

#define BI_COMPLEX_SUMM_THE_PROD(A,B,C);		\
  {							\
    complex_summ_thppe_prod(A[0],B[0],C[0]);		\
    complex_summ_thppe_prod(A[1],B[1],C[1]);		\
  }

#define BI_SU3_DAG_PROD_BI_COLOR(OUT,U,IN)			\
  {								\
    CHECK_BI_SU3_DIFF_FROM_ZERO(U);				\
    BI_COMPLEX_CONJ1_PROD(OUT[0],U[0][0],IN[0]);		\
    BI_COMPLEX_CONJ1_PROD(OUT[1],U[0][1],IN[0]);		\
    BI_COMPLEX_CONJ1_PROD(OUT[2],U[0][2],IN[0]);		\
    BI_COMPLEX_SUMM_THE_CONJ1_PROD(OUT[0],U[1][0],IN[1]);	\
    BI_COMPLEX_SUMM_THE_CONJ1_PROD(OUT[1],U[1][1],IN[1]);	\
    BI_COMPLEX_SUMM_THE_CONJ1_PROD(OUT[2],U[1][2],IN[1]);	\
    BI_COMPLEX_SUMM_THE_CONJ1_PROD(OUT[0],U[2][0],IN[2]);	\
    BI_COMPLEX_SUMM_THE_CONJ1_PROD(OUT[1],U[2][1],IN[2]);	\
    BI_COMPLEX_SUMM_THE_CONJ1_PROD(OUT[2],U[2][2],IN[2]);	\
  }
    
#define BI_SU3_PROD_BI_COLOR(OUT,U,IN)			\
  {							\
    CHECK_BI_SU3_DIFF_FROM_ZERO(U);			\
    BI_COMPLEX_PROD(OUT[0],U[0][0],IN[0]);		\
    BI_COMPLEX_PROD(OUT[1],U[1][0],IN[0]);		\
    BI_COMPLEX_PROD(OUT[2],U[2][0],IN[0]);		\
    BI_COMPLEX_SUMM_THE_PROD(OUT[0],U[0][1],IN[1]);	\
    BI_COMPLEX_SUMM_THE_PROD(OUT[1],U[1][1],IN[1]);	\
    BI_COMPLEX_SUMM_THE_PROD(OUT[2],U[2][1],IN[1]);	\
    BI_COMPLEX_SUMM_THE_PROD(OUT[0],U[0][2],IN[2]);	\
    BI_COMPLEX_SUMM_THE_PROD(OUT[1],U[1][2],IN[2]);	\
    BI_COMPLEX_SUMM_THE_PROD(OUT[2],U[2][2],IN[2]);	\
  }
    
#define BI_SU3_DAG_PROD_BI_HALFSPINCOLOR(A,B,C)	\
  {						\
    BI_SU3_DAG_PROD_BI_COLOR(A[0],B,C[0]);	\
    BI_SU3_DAG_PROD_BI_COLOR(A[1],B,C[1]);	\
  }

#define BI_SU3_PROD_BI_HALFSPINCOLOR(A,B,C)	\
  {						\
    BI_SU3_PROD_BI_COLOR(A[0],B,C[0]);	\
    BI_SU3_PROD_BI_COLOR(A[1],B,C[1]);	\
  }

///////////////////////////////////////////// projectors /////////////////////////////////////////

#define HOPMATR_TBW_PROJ(OUT,IN)		\
  {						\
    BI_COLOR_SUBT(OUT[0],IN[0],IN[2]);		\
    BI_COLOR_SUBT(OUT[1],IN[1],IN[3]);		\
  }

#define HOPMATR_XBW_PROJ(OUT,IN)		\
  {						\
    BI_COLOR_ISUBT(OUT[0],IN[0],IN[3]);		\
    BI_COLOR_ISUBT(OUT[1],IN[1],IN[2]);		\
  }

#define HOPMATR_YBW_PROJ(OUT,IN)		\
  {						\
    BI_COLOR_SUBT(OUT[0],IN[0],IN[3]);		\
    BI_COLOR_SUMM(OUT[1],IN[1],IN[2]);		\
  }

#define HOPMATR_ZBW_PROJ(OUT,IN)		\
  {						\
    BI_COLOR_ISUBT(OUT[0],IN[0],IN[2]);		\
    BI_COLOR_ISUMM(OUT[1],IN[1],IN[3]);		\
  }

#define HOPMATR_TFW_PROJ(OUT,IN)		\
  {						\
    BI_COLOR_SUMM(OUT[0],IN[0],IN[2]);		\
    BI_COLOR_SUMM(OUT[1],IN[1],IN[3]);		\
  }

#define HOPMATR_XFW_PROJ(OUT,IN)		\
  {						\
    BI_COLOR_ISUMM(OUT[0],IN[0],IN[3]);		\
    BI_COLOR_ISUMM(OUT[1],IN[1],IN[2]);		\
  }

#define HOPMATR_YFW_PROJ(OUT,IN)		\
  {						\
    BI_COLOR_SUMM(OUT[0],IN[0],IN[3]);		\
    BI_COLOR_SUBT(OUT[1],IN[1],IN[2]);		\
  }

#define HOPMATR_ZFW_PROJ(OUT,IN)		\
  {						\
    BI_COLOR_ISUMM(OUT[0],IN[0],IN[2]);		\
    BI_COLOR_ISUBT(OUT[1],IN[1],IN[3]);		\
  }

/////////////////////////////////////////////////////////////////////////////////////

#define BI_COMPLEX_TRANSPOSE(A,B)		\
  {						\
    A[0][0]=B[1][0];				\
    A[0][1]=B[1][1];				\
    A[1][0]=B[0][0];				\
    A[1][1]=B[0][1];				\
  }

#define BI_COLOR_TRANSPOSE(A,B)			\
  {						\
    BI_COMPLEX_TRANSPOSE(A[0],B[0]);		\
    BI_COMPLEX_TRANSPOSE(A[1],B[1]);		\
    BI_COMPLEX_TRANSPOSE(A[2],B[2]);		\
  }

#define BI_HALFSPINCOLOR_TRANSPOSE(A,B)		\
  {						\
    BI_COLOR_TRANSPOSE(A[0],B[0]);		\
    BI_COLOR_TRANSPOSE(A[1],B[1]);		\
  }

///////////////////////////////////////////// expand hppalfspincolors ///////////////////////////////////////

#define DIAG_TMQ(OUT,DIAG,IN)				\
  {							\
    BI_COLOR_PROD_COMPLEX(OUT[0],IN[0],DIAG[0]);	\
    BI_COLOR_PROD_COMPLEX(OUT[1],IN[1],DIAG[0]);	\
    BI_COLOR_PROD_COMPLEX(OUT[2],IN[2],DIAG[1]);	\
    BI_COLOR_PROD_COMPLEX(OUT[3],IN[3],DIAG[1]);	\
  }

#define TFW_DER_TMQ_EXP(OUT,IN)				\
  {							\
    BI_COLOR_SUMMASSIGN(OUT[0],IN[0]);			\
    BI_COLOR_SUMMASSIGN(OUT[1],IN[1]);			\
    BI_COLOR_SUBTASSIGN(OUT[2],IN[0]);			\
    BI_COLOR_SUBTASSIGN(OUT[3],IN[1]);			\
  }

#define XFW_DER_TMQ_EXP(OUT,IN)				\
  {							\
    BI_COLOR_SUMMASSIGN(OUT[0],IN[0]);			\
    BI_COLOR_SUMMASSIGN(OUT[1],IN[1]);			\
    BI_COLOR_ISUMMASSIGN(OUT[2],IN[1]);			\
    BI_COLOR_ISUMMASSIGN(OUT[3],IN[0]);			\
  }

#define YFW_DER_TMQ_EXP(OUT,IN)				\
  {							\
    BI_COLOR_SUMMASSIGN(OUT[0],IN[0]);			\
    BI_COLOR_SUMMASSIGN(OUT[1],IN[1]);			\
    BI_COLOR_SUMMASSIGN(OUT[2],IN[1]);			\
    BI_COLOR_SUBTASSIGN(OUT[3],IN[0]);			\
  }

#define ZFW_DER_TMQ_EXP(OUT,IN)				\
  {							\
    BI_COLOR_SUMMASSIGN(OUT[0],IN[0]);			\
    BI_COLOR_SUMMASSIGN(OUT[1],IN[1]);			\
    BI_COLOR_ISUMMASSIGN(OUT[2],IN[0]);			\
    BI_COLOR_ISUBTASSIGN(OUT[3],IN[1]);			\
  }

#define TBW_DER_TMQ_EXP(OUT,IN)				\
  {							\
    BI_COLOR_SUMMASSIGN(OUT[0],IN[0]);			\
    BI_COLOR_SUMMASSIGN(OUT[1],IN[1]);			\
    BI_COLOR_SUMMASSIGN(OUT[2],IN[0]);			\
    BI_COLOR_SUMMASSIGN(OUT[3],IN[1]);			\
  }

#define XBW_DER_TMQ_EXP(OUT,IN)				\
  {							\
    BI_COLOR_SUMMASSIGN(OUT[0],IN[0]);			\
    BI_COLOR_SUMMASSIGN(OUT[1],IN[1]);			\
    BI_COLOR_ISUBTASSIGN(OUT[2],IN[1]);			\
    BI_COLOR_ISUBTASSIGN(OUT[3],IN[0]);			\
  }

#define YBW_DER_TMQ_EXP(OUT,IN)				\
  {							\
    BI_COLOR_SUMMASSIGN(OUT[0],IN[0]);			\
    BI_COLOR_SUMMASSIGN(OUT[1],IN[1]);			\
    BI_COLOR_SUBTASSIGN(OUT[2],IN[1]);			\
    BI_COLOR_SUMMASSIGN(OUT[3],IN[0]);			\
  }

#define ZBW_DER_TMQ_EXP(OUT,IN)				\
  {							\
    BI_COLOR_SUMMASSIGN(OUT[0],IN[0]);			\
    BI_COLOR_SUMMASSIGN(OUT[1],IN[1]);			\
    BI_COLOR_ISUBTASSIGN(OUT[2],IN[0]);			\
    BI_COLOR_ISUMMASSIGN(OUT[3],IN[1]);			\
  }

