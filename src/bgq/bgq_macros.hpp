#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "intrinsic.hpp"
#include "new_types/complex.hpp"
#include "new_types/float_128.hpp"

#define BI_SU3_PUT_TO_ZERO(A) memset(A,0,sizeof(bi_su3))

////////////////////////////// convert normal complex to BI //////////////////////////

#ifdef BGQ_EMU
 #define DECLARE_REG_BI_COMPLEX(A) bi_complex A
 #define COMPLEX_TO_BI_COMPLEX(A,B,VN) complex_copy(A[VN],B)
 #define COMPLEX_TO_BI_SINGLE_COMPLEX(A,B,VN) single_complex_copy_from_complex(A[VN],B)
#else
#define DECLARE_REG_BI_COMPLEX(A) reg_bi_complex A
 #define COMPLEX_TO_BI_COMPLEX(A,B,VN) vec_st2(vec_ld2(0,B),0,A[VN])
 #define COMPLEX_TO_BI_SINGLE_COMPLEX(A,B,VN) vec_st2(vec_ld2(0,B),0,A[VN])
#endif

//maybe we will promote it
#define COMPLEX_TO_BI_COMPLEX_CONJ(A,B,VN) complex_conj(A[VN],B)

#define COLOR_TO_BI_COLOR(A,B,VN)		\
  do						\
    {						\
      COMPLEX_TO_BI_COMPLEX(A[0],B[0],VN);	\
      COMPLEX_TO_BI_COMPLEX(A[1],B[1],VN);	\
      COMPLEX_TO_BI_COMPLEX(A[2],B[2],VN);	\
    }						\
  while(0)
#define COLOR_TO_BI_SINGLE_COLOR(A,B,VN)	\
  do							\
    {							\
      COMPLEX_TO_BI_SINGLE_COMPLEX(A[0],B[0],VN);	\
      COMPLEX_TO_BI_SINGLE_COMPLEX(A[1],B[1],VN);	\
      COMPLEX_TO_BI_SINGLE_COMPLEX(A[2],B[2],VN);	\
    }							\
  while(0)
#define SU3_TO_BI_SU3(A,B,VN)			\
  do						\
    {						\
      COLOR_TO_BI_COLOR(A[0],B[0],VN);		\
      COLOR_TO_BI_COLOR(A[1],B[1],VN);		\
      COLOR_TO_BI_COLOR(A[2],B[2],VN);		\
    }						\
  while(0)
#define SU3_TO_BI_SU3_DAG(A,B,VN)					\
  for(int icol=0;icol<NCOL;icol++)					\
    for(int jcol=0;jcol<NCOL;jcol++)					\
      COMPLEX_TO_BI_COMPLEX_CONJ(A[jcol][icol],B[icol][jcol],VN)
#define SU3_TO_BI_SINGLE_SU3(A,B,VN)		\
  do						\
    {						\
      COLOR_TO_BI_SINGLE_COLOR(A[0],B[0],VN);	\
      COLOR_TO_BI_SINGLE_COLOR(A[1],B[1],VN);	\
      COLOR_TO_BI_SINGLE_COLOR(A[2],B[2],VN);	\
    }							\
  while(0)
#define HALFSPINCOLOR_TO_BI_HALFSPINCOLOR(A,B,VN)	\
  do							\
    {							\
      COLOR_TO_BI_COLOR(A[0],B[0],VN);			\
      COLOR_TO_BI_COLOR(A[1],B[1],VN);			\
    }							\
  while(0)
#define SPINCOLOR_TO_BI_SPINCOLOR(A,B,VN)	\
  do						\
    {						\
      COLOR_TO_BI_COLOR(A[0],B[0],VN);		\
      COLOR_TO_BI_COLOR(A[1],B[1],VN);		\
      COLOR_TO_BI_COLOR(A[2],B[2],VN);		\
      COLOR_TO_BI_COLOR(A[3],B[3],VN);		\
    }						\
  while(0)

/////////////// 128 bit version //////////////

#define COMPLEX_128_TO_BI_COMPLEX_128(A,B,VN)		\
  do							\
    {							\
      float_128_copy(A[VN][RE],B[RE]);			\
      float_128_copy(A[VN][IM],B[IM]);			\
    }							\
  while(0)
#define COLOR_128_TO_BI_COLOR_128(A,B,VN)		\
  do							\
    {							\
      COMPLEX_128_TO_BI_COMPLEX_128(A[0],B[0],VN);	\
      COMPLEX_128_TO_BI_COMPLEX_128(A[1],B[1],VN);	\
      COMPLEX_128_TO_BI_COMPLEX_128(A[2],B[2],VN);	\
    }							\
  while(0)
#define SU3_128_TO_BI_SU3_128(A,B,VN)		\
  do						\
    {						\
      COLOR_128_TO_BI_COLOR_128(A[0],B[0],VN);	\
      COLOR_128_TO_BI_COLOR_128(A[1],B[1],VN);	\
      COLOR_128_TO_BI_COLOR_128(A[2],B[2],VN);	\
    }								\
  while(0)
#define HALFSPINCOLOR_128_TO_BI_HALFSPINCOLOR_128(A,B,VN)	\
  do								\
    {								\
      COLOR_128_TO_BI_COLOR_128(A[0],B[0],VN);			\
      COLOR_128_TO_BI_COLOR_128(A[1],B[1],VN);			\
    }								\
  while(0)
#define SPINCOLOR_128_TO_BI_SPINCOLOR_128(A,B,VN)	\
  do							\
    {							\
      COLOR_128_TO_BI_COLOR_128(A[0],B[0],VN);		\
      COLOR_128_TO_BI_COLOR_128(A[1],B[1],VN);		\
      COLOR_128_TO_BI_COLOR_128(A[2],B[2],VN);		\
      COLOR_128_TO_BI_COLOR_128(A[3],B[3],VN);		\
    }							\
  while(0)

//////////////////////////////// extract component of BI ////////////////////////////////

#define COMPLEX_OF_BI_COMPLEX(A,B,VN)		\
    complex_copy(A,B[VN])
#define COLOR_OF_BI_COLOR(A,B,VN)		\
  do						\
    {						\
      COMPLEX_OF_BI_COMPLEX(A[0],B[0],VN);	\
      COMPLEX_OF_BI_COMPLEX(A[1],B[1],VN);	\
      COMPLEX_OF_BI_COMPLEX(A[2],B[2],VN);	\
    }						\
  while(0)
#define SU3_OF_BI_SU3(A,B,VN)			\
  do						\
    {						\
      COLOR_OF_BI_COLOR(A[0],B[0],VN);		\
      COLOR_OF_BI_COLOR(A[1],B[1],VN);		\
      COLOR_OF_BI_COLOR(A[2],B[2],VN);		\
    }							\
  while(0)
#define HALFSPINCOLOR_OF_BI_HALFSPINCOLOR(A,B,VN)	\
  do							\
    {							\
      COLOR_OF_BI_COLOR(A[0],B[0],VN);			\
      COLOR_OF_BI_COLOR(A[1],B[1],VN);			\
    }							\
  while(0)
#define SPINCOLOR_OF_BI_SPINCOLOR(A,B,VN)	\
  do						\
    {						\
      COLOR_OF_BI_COLOR(A[0],B[0],VN);		\
      COLOR_OF_BI_COLOR(A[1],B[1],VN);		\
      COLOR_OF_BI_COLOR(A[2],B[2],VN);		\
      COLOR_OF_BI_COLOR(A[3],B[3],VN);		\
    }						\
  while(0)

///////////// 128 bit version /////////////

#define COMPLEX_128_OF_BI_COMPLEX_128(A,B,VN)		\
  do							\
    {							\
      float_128_copy(A[RE],B[VN][RE]);			\
      float_128_copy(A[IM],B[VN][IM]);			\
    }							\
  while(0)
#define COLOR_128_OF_BI_COLOR_128(A,B,VN)		\
  do							\
    {							\
      COMPLEX_128_OF_BI_COMPLEX_128(A[0],B[0],VN);	\
      COMPLEX_128_OF_BI_COMPLEX_128(A[1],B[1],VN);	\
      COMPLEX_128_OF_BI_COMPLEX_128(A[2],B[2],VN);	\
    }							\
  while(0)
#define SU3_128_OF_BI_SU3_128(A,B,VN)		\
  do						\
    {						\
      COLOR_128_OF_BI_COLOR_128(A[0],B[0],VN);	\
      COLOR_128_OF_BI_COLOR_128(A[1],B[1],VN);	\
      COLOR_128_OF_BI_COLOR_128(A[2],B[2],VN);	\
    }								\
  while(0)
#define HALFSPINCOLOR_128_OF_BI_HALFSPINCOLOR_128(A,B,VN)	\
  do								\
    {								\
      COLOR_128_OF_BI_COLOR_128(A[0],B[0],VN);			\
      COLOR_128_OF_BI_COLOR_128(A[1],B[1],VN);			\
    }								\
  while(0)
#define SPINCOLOR_128_OF_BI_SPINCOLOR_128(A,B,VN)	\
  do							\
    {							\
      COLOR_128_OF_BI_COLOR_128(A[0],B[0],VN);		\
      COLOR_128_OF_BI_COLOR_128(A[1],B[1],VN);		\
      COLOR_128_OF_BI_COLOR_128(A[2],B[2],VN);		\
      COLOR_128_OF_BI_COLOR_128(A[3],B[3],VN);		\
    }							\
  while(0)

//////////////////////////////// copy BI ////////////////////////////////

#define BI_COMPLEX_SPLAT(A,B) A[0][0]=A[0][1]=A[1][0]=A[1][1]=B

#define BI_COMPLEX_COPY(A,B)			\
  do						\
    {						\
    complex_copy(A[0],B[0]);			\
    complex_copy(A[1],B[1]);			\
    }						\
  while(0)

#define BI_COMPLEX_COPY_FROM_BI_SINGLE_COMPLEX(A,B)	\
  do							\
  {							\
    complex_copy_from_single_complex(A[0],B[0]);	\
    complex_copy_from_single_complex(A[1],B[1]);	\
  }							\
  while(0)
#define BI_SINGLE_COMPLEX_COPY_FROM_BI_COMPLEX(A,B)	\
  do							\
  {							\
    single_complex_copy_from_complex(A[0],B[0]);	\
    single_complex_copy_from_complex(A[1],B[1]);	\
  }							\
  while(0)
#define BI_COLOR_COPY(A,B)			\
  do						\
  {						\
    BI_COMPLEX_COPY(A[0],B[0]);			\
    BI_COMPLEX_COPY(A[1],B[1]);			\
    BI_COMPLEX_COPY(A[2],B[2]);			\
  }						\
  while(0)
#define BI_SU3_COPY(A,B)			\
  do						\
  {						\
    BI_COLOR_COPY(A[0],B[0]);			\
    BI_COLOR_COPY(A[1],B[1]);			\
    BI_COLOR_COPY(A[2],B[2]);			\
  }						\
  while(0)
#define BI_HALFSPINCOLOR_COPY(A,B)		\
  do						\
  {						\
    BI_COLOR_COPY(A[0],B[0]);			\
    BI_COLOR_COPY(A[1],B[1]);			\
  }						\
  while(0)
#define BI_SPINCOLOR_COPY(A,B)			\
  do						\
  {						\
    BI_COLOR_COPY(A[0],B[0]);			\
    BI_COLOR_COPY(A[1],B[1]);			\
    BI_COLOR_COPY(A[2],B[2]);			\
    BI_COLOR_COPY(A[3],B[3]);			\
  }						\
  while(0)

/////////////////////////////////// split BI ///////////////////////////////////

#define BI_COMPLEX_TO_COMPLEX(A,B,C)		\
  do						\
    {						\
      complex_copy(A,C[0]);			\
      complex_copy(B,C[1]);			\
    }						\
  while(0)
#define BI_SINGLE_COMPLEX_TO_COMPLEX(A,B,C)	\
  do						\
    {						\
      complex_copy_from_single_complex(A,C[0]);	\
      complex_copy_from_single_complex(B,C[1]);	\
    }						\
  while(0)
#define BI_COLOR_TO_COLOR(A,B,C)		\
  do						\
    {						\
      BI_COMPLEX_TO_COMPLEX(A[0],B[0],C[0]);	\
      BI_COMPLEX_TO_COMPLEX(A[1],B[1],C[1]);	\
      BI_COMPLEX_TO_COMPLEX(A[2],B[2],C[2]);	\
    }							\
  while(0)
#define BI_SINGLE_COLOR_TO_COLOR(A,B,C)			\
  do							\
    {							\
      BI_SINGLE_COMPLEX_TO_COMPLEX(A[0],B[0],C[0]);	\
      BI_SINGLE_COMPLEX_TO_COMPLEX(A[1],B[1],C[1]);	\
      BI_SINGLE_COMPLEX_TO_COMPLEX(A[2],B[2],C[2]);	\
    }							\
  while(0)
#define BI_SU3_TO_SU3(A,B,C)			\
  do						\
    {						\
      BI_COLOR_TO_COLOR(A[0],B[0],C[0]);	\
      BI_COLOR_TO_COLOR(A[1],B[1],C[1]);	\
      BI_COLOR_TO_COLOR(A[2],B[2],C[2]);	\
    }							\
  while(0)
#define BI_HALFSPINCOLOR_TO_HALFSPINCOLOR(A,B,C)	\
  do							\
    {							\
      BI_COLOR_TO_COLOR(A[0],B[0],C[0]);		\
      BI_COLOR_TO_COLOR(A[1],B[1],C[1]);		\
    }							\
  while(0)
#define BI_SPINCOLOR_TO_SPINCOLOR(A,B,C)	\
  do						\
    {						\
      BI_COLOR_TO_COLOR(A[0],B[0],C[0]);	\
      BI_COLOR_TO_COLOR(A[1],B[1],C[1]);	\
      BI_COLOR_TO_COLOR(A[2],B[2],C[2]);	\
      BI_COLOR_TO_COLOR(A[3],B[3],C[3]);	\
    }						\
  while(0)


///////////////// 128 bit version ////////////////

#define BI_COMPLEX_128_TO_COMPLEX_128(A,B,C)		\
  do							\
    {							\
      float_128_copy(A[RE],C[0][RE]);			\
      float_128_copy(A[IM],C[0][IM]);			\
      float_128_copy(B[RE],C[1][RE]);			\
      float_128_copy(B[IM],C[1][IM]);			\
    }							\
  while(0)
#define BI_COLOR_128_TO_COLOR_128(A,B,C)		\
  do							\
    {							\
      BI_COMPLEX_128_TO_COMPLEX_128(A[0],B[0],C[0]);	\
      BI_COMPLEX_128_TO_COMPLEX_128(A[1],B[1],C[1]);	\
      BI_COMPLEX_128_TO_COMPLEX_128(A[2],B[2],C[2]);	\
    }							\
  while(0)
#define BI_SU3_128_TO_SU3_128(A,B,C)		\
  do							\
    {							\
      BI_COLOR_128_TO_COLOR_128(A[0],B[0],C[0]);	\
      BI_COLOR_128_TO_COLOR_128(A[1],B[1],C[1]);	\
      BI_COLOR_128_TO_COLOR_128(A[2],B[2],C[2]);	\
    }								\
  while(0)
#define BI_HALFSPINCOLOR_128_TO_HALFSPINCOLOR_128(A,B,C)	\
  do								\
    {								\
      BI_COLOR_128_TO_COLOR_128(A[0],B[0],C[0]);		\
      BI_COLOR_128_TO_COLOR_128(A[1],B[1],C[1]);		\
    }								\
  while(0)
#define BI_SPINCOLOR_128_TO_SPINCOLOR_128(A,B,C)	\
  do							\
    {							\
      BI_COLOR_128_TO_COLOR_128(A[0],B[0],C[0]);	\
      BI_COLOR_128_TO_COLOR_128(A[1],B[1],C[1]);	\
      BI_COLOR_128_TO_COLOR_128(A[2],B[2],C[2]);	\
      BI_COLOR_128_TO_COLOR_128(A[3],B[3],C[3]);	\
    }							\
  while(0)

//////////////////////////////////// basic operation on BI ////////////////////////////////

#define BI_COMPLEX_PRINT(A)			\
  do						\
    {						\
      complex_print(A[0]);			\
      complex_print(A[1]);			\
    }						\
  while(0)

#define BI_COMPLEX_SUBT(A,B,C)                    \
  do						  \
    {						  \
      complex_subt(A[0],B[0],C[0]);		  \
      complex_subt(A[1],B[1],C[1]);		  \
    }						  \
  while(0)
#define BI_COMPLEX_SUMM(A,B,C)                    \
  do						  \
    {						  \
      complex_summ(A[0],B[0],C[0]);		  \
      complex_summ(A[1],B[1],C[1]);		  \
    }						  \
  while(0)

#define BI_COMPLEX_ISUBT(A,B,C)			  \
  do						  \
    {						  \
      complex_isubt(A[0],B[0],C[0]);		  \
      complex_isubt(A[1],B[1],C[1]);		  \
    }						  \
  while(0)

#define BI_COMPLEX_ISUMM(A,B,C)			  \
  do						  \
    {						  \
      complex_isumm(A[0],B[0],C[0]);		  \
      complex_isumm(A[1],B[1],C[1]);		  \
    }						  \
  while(0)

#define BI_COLOR_SUBT(A,B,C)			\
  do						\
    {						\
      BI_COMPLEX_SUBT(A[0],B[0],C[0]);		\
      BI_COMPLEX_SUBT(A[1],B[1],C[1]);		\
      BI_COMPLEX_SUBT(A[2],B[2],C[2]);		\
    }						\
  while(0)

#define BI_COLOR_SUMM(A,B,C)			\
  do						\
    {						\
      BI_COMPLEX_SUMM(A[0],B[0],C[0]);		\
      BI_COMPLEX_SUMM(A[1],B[1],C[1]);		\
      BI_COMPLEX_SUMM(A[2],B[2],C[2]);		\
    }						\
  while(0)

#define BI_COLOR_ISUBT(A,B,C)			\
  do						\
    {						\
      BI_COMPLEX_ISUBT(A[0],B[0],C[0]);		\
      BI_COMPLEX_ISUBT(A[1],B[1],C[1]);		\
      BI_COMPLEX_ISUBT(A[2],B[2],C[2]);		\
    }						\
  while(0)

#define BI_COLOR_ISUMM(A,B,C)			\
  do						\
    {						\
      BI_COMPLEX_ISUMM(A[0],B[0],C[0]);		\
      BI_COMPLEX_ISUMM(A[1],B[1],C[1]);		\
      BI_COMPLEX_ISUMM(A[2],B[2],C[2]);		\
    }						\
  while(0)

#define BI_COLOR_PROD_COMPLEX(A,B,C)		\
  do						\
    {						\
      BI_COMPLEX_PROD(A[0],B[0],C);		\
      BI_COMPLEX_PROD(A[1],B[1],C);		\
      BI_COMPLEX_PROD(A[2],B[2],C);		\
    }						\
  while(0)

#define BI_COLOR_PROD_DOUBLE(A,B,C)		\
  do						\
    {						\
      BI_COMPLEX_PROD_DOUBLE(A[0],B[0],C);	\
      BI_COMPLEX_PROD_DOUBLE(A[1],B[1],C);	\
      BI_COMPLEX_PROD_DOUBLE(A[2],B[2],C);	\
    }						\
  while(0)

#define BI_COMPLEX_PROD_4DOUBLE(A,B,C)		\
  do						\
    {						\
      A[0][0]=B[0][0]*C[0][0];			\
      A[0][1]=B[0][1]*C[0][1];			\
      A[1][0]=B[1][0]*C[1][0];			\
      A[1][1]=B[1][1]*C[1][1];			\
    }						\
  while(0)

#define BI_COLOR_SUMMASSIGN(A,B) BI_COLOR_SUMM(A,A,B)
#define BI_COLOR_SUBTASSIGN(A,B) BI_COLOR_SUBT(A,A,B)
#define BI_COLOR_ISUMMASSIGN(A,B) BI_COLOR_ISUMM(A,A,B)
#define BI_COLOR_ISUBTASSIGN(A,B) BI_COLOR_ISUBT(A,A,B)

//////////////////////////////////////////////////////////////////////////////////////////////////////

#define BI_HALFSPINCOLOR_SUMM(A,B,C)		\
  do						\
    {						\
      BI_COLOR_SUMM(A[0],B[0],C[0]);		\
      BI_COLOR_SUMM(A[1],B[1],C[1]);		\
    }						\
  while(0)

#define BI_HALFSPINCOLOR_SUMMASSIGN(A,B)	\
  BI_HALFSPINCOLOR_SUMM(A,A,B);

#define BI_SPINCOLOR_PROD_DOUBLE(A,B,C)					\
  do									\
    {									\
      BI_COLOR_PROD_DOUBLE(A[0],B[0],C);				\
      BI_COLOR_PROD_DOUBLE(A[1],B[1],C);				\
      BI_COLOR_PROD_DOUBLE(A[2],B[2],C);				\
      BI_COLOR_PROD_DOUBLE(A[3],B[3],C);				\
    }									\
  while(0)

//////////////////////////////////////////////////////////////////////////////////////////////////////

#define BI_COMPLEX_CONJ1_PROD(A,B,C)			\
  do							\
    {							\
      unsafe_complex_conj1_prod(A[0],B[0],C[0]);	\
      unsafe_complex_conj1_prod(A[1],B[1],C[1]);	\
    }							\
  while(0)

#define BI_COMPLEX_PROD(A,B,C)			\
  do						\
    {						\
      unsafe_complex_prod(A[0],B[0],C[0]);	\
      unsafe_complex_prod(A[1],B[1],C[1]);	\
    }						\
  while(0)

#define BI_COMPLEX_PROD_DOUBLE(A,B,C)		\
  do						\
    {						\
      complex_prod_double(A[0],B[0],C);		\
      complex_prod_double(A[1],B[1],C);		\
    }						\
  while(0)

#define BI_COMPLEX_SUMM_THE_CONJ1_PROD(A,B,C)		\
  do							\
    {							\
      complex_summ_the_conj1_prod(A[0],B[0],C[0]);	\
      complex_summ_the_conj1_prod(A[1],B[1],C[1]);	\
    }							\
  while(0)

#define BI_COMPLEX_SUBT_THE_CONJ1_PROD(A,B,C)		\
  do							\
    {							\
      complex_subt_the_conj1_prod(A[0],B[0],C[0]);	\
      complex_subt_the_conj1_prod(A[1],B[1],C[1]);	\
    }							\
  while(0)

#define BI_COMPLEX_SUMM_THE_PROD(A,B,C)			\
  do							\
    {							\
      complex_summ_the_prod(A[0],B[0],C[0]);		\
      complex_summ_the_prod(A[1],B[1],C[1]);		\
    }							\
  while(0)

#define BI_COMPLEX_SUBT_THE_PROD(A,B,C)			\
  do							\
    {							\
      complex_subt_the_prod(A[0],B[0],C[0]);		\
      complex_subt_the_prod(A[1],B[1],C[1]);		\
    }							\
  while(0)

#define BI_COMPLEX_SUMM_THE_PROD_4DOUBLE(A,B,C,D)	\
  do							\
    {							\
      A[0][0]=B[0][0]+C[0][0]*D[0][0];			\
      A[0][1]=B[0][1]+C[0][1]*D[0][1];			\
      A[1][0]=B[1][0]+C[1][0]*D[1][0];			\
      A[1][1]=B[1][1]+C[1][1]*D[1][1];			\
    }							\
  while(0)

#define BI_SU3_DAG_PROD_BI_COLOR(OUT,U,IN)			\
  do								\
    {								\
      BI_COMPLEX_CONJ1_PROD(OUT[0],U[0][0],IN[0]);		\
      BI_COMPLEX_CONJ1_PROD(OUT[1],U[0][1],IN[0]);		\
      BI_COMPLEX_CONJ1_PROD(OUT[2],U[0][2],IN[0]);		\
      BI_COMPLEX_SUMM_THE_CONJ1_PROD(OUT[0],U[1][0],IN[1]);	\
      BI_COMPLEX_SUMM_THE_CONJ1_PROD(OUT[1],U[1][1],IN[1]);	\
      BI_COMPLEX_SUMM_THE_CONJ1_PROD(OUT[2],U[1][2],IN[1]);	\
      BI_COMPLEX_SUMM_THE_CONJ1_PROD(OUT[0],U[2][0],IN[2]);	\
      BI_COMPLEX_SUMM_THE_CONJ1_PROD(OUT[1],U[2][1],IN[2]);	\
      BI_COMPLEX_SUMM_THE_CONJ1_PROD(OUT[2],U[2][2],IN[2]);	\
    }								\
  while(0)

#define BI_SU3_PROD_BI_COLOR(OUT,U,IN)			\
  do							\
    {							\
      BI_COMPLEX_PROD(OUT[0],U[0][0],IN[0]);		\
      BI_COMPLEX_PROD(OUT[1],U[1][0],IN[0]);		\
      BI_COMPLEX_PROD(OUT[2],U[2][0],IN[0]);		\
      BI_COMPLEX_SUMM_THE_PROD(OUT[0],U[0][1],IN[1]);	\
      BI_COMPLEX_SUMM_THE_PROD(OUT[1],U[1][1],IN[1]);	\
      BI_COMPLEX_SUMM_THE_PROD(OUT[2],U[2][1],IN[1]);	\
      BI_COMPLEX_SUMM_THE_PROD(OUT[0],U[0][2],IN[2]);	\
      BI_COMPLEX_SUMM_THE_PROD(OUT[1],U[1][2],IN[2]);	\
      BI_COMPLEX_SUMM_THE_PROD(OUT[2],U[2][2],IN[2]);	\
    }							\
  while(0)

#define BI_SU3_DAG_PROD_BI_HALFSPINCOLOR(A,B,C)	\
  do						\
    {						\
      BI_SU3_DAG_PROD_BI_COLOR(A[0],B,C[0]);	\
      BI_SU3_DAG_PROD_BI_COLOR(A[1],B,C[1]);	\
    }

#define BI_SU3_PROD_BI_HALFSPINCOLOR(A,B,C)	\
  do						\
    {						\
      BI_SU3_PROD_BI_COLOR(A[0],B,C[0]);	\
      BI_SU3_PROD_BI_COLOR(A[1],B,C[1]);	\
    }

///////////////////////////////////////////// projectors /////////////////////////////////////////

#define HOPMATR_TBW_PROJ(OUT,IN)		\
  do						\
    {						\
      BI_COLOR_SUBT(OUT[0],IN[0],IN[2]);	\
      BI_COLOR_SUBT(OUT[1],IN[1],IN[3]);	\
    }

#define HOPMATR_XBW_PROJ(OUT,IN)		\
  do						\
    {						\
      BI_COLOR_ISUBT(OUT[0],IN[0],IN[3]);	\
      BI_COLOR_ISUBT(OUT[1],IN[1],IN[2]);	\
    }						\
  while(0)

#define HOPMATR_YBW_PROJ(OUT,IN)		\
  do						\
    {						\
      BI_COLOR_SUBT(OUT[0],IN[0],IN[3]);	\
      BI_COLOR_SUMM(OUT[1],IN[1],IN[2]);	\
    }						\
  while(0)

#define HOPMATR_ZBW_PROJ(OUT,IN)		\
  do						\
    {						\
      BI_COLOR_ISUBT(OUT[0],IN[0],IN[2]);	\
      BI_COLOR_ISUMM(OUT[1],IN[1],IN[3]);	\
    }						\
  while(0)

#define HOPMATR_TFW_PROJ(OUT,IN)		\
  do						\
    {						\
      BI_COLOR_SUMM(OUT[0],IN[0],IN[2]);	\
      BI_COLOR_SUMM(OUT[1],IN[1],IN[3]);	\
    }						\
  while(0)

#define HOPMATR_XFW_PROJ(OUT,IN)		\
  do						\
    {						\
      BI_COLOR_ISUMM(OUT[0],IN[0],IN[3]);	\
      BI_COLOR_ISUMM(OUT[1],IN[1],IN[2]);	\
    }						\
  while(0)

#define HOPMATR_YFW_PROJ(OUT,IN)		\
  do						\
    {						\
      BI_COLOR_SUMM(OUT[0],IN[0],IN[3]);	\
      BI_COLOR_SUBT(OUT[1],IN[1],IN[2]);	\
    }						\
  while(0)

#define HOPMATR_ZFW_PROJ(OUT,IN)		\
  do						\
    {						\
      BI_COLOR_ISUMM(OUT[0],IN[0],IN[2]);	\
      BI_COLOR_ISUBT(OUT[1],IN[1],IN[3]);	\
    }						\
  while(0)

/////////////////////////////////////////////////////////////////////////////////////

#define BI_COMPLEX_TRANSPOSE(A,B)		\
  do						\
    {						\
      A[0][0]=B[1][0];				\
      A[0][1]=B[1][1];				\
      A[1][0]=B[0][0];				\
      A[1][1]=B[0][1];				\
    }						\
  while(0)

#define BI_COLOR_TRANSPOSE(A,B)			\
  do						\
    {						\
      BI_COMPLEX_TRANSPOSE(A[0],B[0]);		\
      BI_COMPLEX_TRANSPOSE(A[1],B[1]);		\
      BI_COMPLEX_TRANSPOSE(A[2],B[2]);		\
    }						\
  while(0)

#define BI_HALFSPINCOLOR_TRANSPOSE(A,B)		\
  do						\
    {						\
      BI_COLOR_TRANSPOSE(A[0],B[0]);		\
      BI_COLOR_TRANSPOSE(A[1],B[1]);		\
    }						\
  while(0)

///////////////////////////////////////////// expand halfspincolors ///////////////////////////////////////

#define DIAG_TMQ(OUT,DIAG,IN)				\
  do							\
    {							\
      BI_COLOR_PROD_COMPLEX(OUT[0],IN[0],DIAG[0]);	\
      BI_COLOR_PROD_COMPLEX(OUT[1],IN[1],DIAG[0]);	\
      BI_COLOR_PROD_COMPLEX(OUT[2],IN[2],DIAG[1]);	\
      BI_COLOR_PROD_COMPLEX(OUT[3],IN[3],DIAG[1]);	\
    }							\
  while(0)

#define TFW_DER_TMQ_EXP(OUT,IN)				\
  do							\
    {							\
      BI_COLOR_SUMMASSIGN(OUT[0],IN[0]);		\
      BI_COLOR_SUMMASSIGN(OUT[1],IN[1]);		\
      BI_COLOR_SUBTASSIGN(OUT[2],IN[0]);		\
      BI_COLOR_SUBTASSIGN(OUT[3],IN[1]);		\
    }							\
  while(0)

#define XFW_DER_TMQ_EXP(OUT,IN)				\
  do							\
    {							\
      BI_COLOR_SUMMASSIGN(OUT[0],IN[0]);		\
      BI_COLOR_SUMMASSIGN(OUT[1],IN[1]);		\
      BI_COLOR_ISUMMASSIGN(OUT[2],IN[1]);		\
      BI_COLOR_ISUMMASSIGN(OUT[3],IN[0]);		\
    }							\
  while(0)

#define YFW_DER_TMQ_EXP(OUT,IN)				\
  do							\
    {							\
      BI_COLOR_SUMMASSIGN(OUT[0],IN[0]);		\
      BI_COLOR_SUMMASSIGN(OUT[1],IN[1]);		\
      BI_COLOR_SUMMASSIGN(OUT[2],IN[1]);		\
      BI_COLOR_SUBTASSIGN(OUT[3],IN[0]);		\
    }							\
  while(0)

#define ZFW_DER_TMQ_EXP(OUT,IN)				\
  do							\
    {							\
      BI_COLOR_SUMMASSIGN(OUT[0],IN[0]);		\
      BI_COLOR_SUMMASSIGN(OUT[1],IN[1]);		\
      BI_COLOR_ISUMMASSIGN(OUT[2],IN[0]);		\
      BI_COLOR_ISUBTASSIGN(OUT[3],IN[1]);		\
    }							\
  while(0)

#define TBW_DER_TMQ_EXP(OUT,IN)				\
  do							\
    {							\
      BI_COLOR_SUMMASSIGN(OUT[0],IN[0]);		\
      BI_COLOR_SUMMASSIGN(OUT[1],IN[1]);		\
      BI_COLOR_SUMMASSIGN(OUT[2],IN[0]);		\
      BI_COLOR_SUMMASSIGN(OUT[3],IN[1]);		\
    }							\
  while(0)

#define XBW_DER_TMQ_EXP(OUT,IN)				\
  do							\
    {							\
      BI_COLOR_SUMMASSIGN(OUT[0],IN[0]);		\
      BI_COLOR_SUMMASSIGN(OUT[1],IN[1]);		\
      BI_COLOR_ISUBTASSIGN(OUT[2],IN[1]);		\
      BI_COLOR_ISUBTASSIGN(OUT[3],IN[0]);		\
    }							\
  while(0)

#define YBW_DER_TMQ_EXP(OUT,IN)				\
  do							\
    {							\
      BI_COLOR_SUMMASSIGN(OUT[0],IN[0]);		\
      BI_COLOR_SUMMASSIGN(OUT[1],IN[1]);		\
      BI_COLOR_SUBTASSIGN(OUT[2],IN[1]);		\
      BI_COLOR_SUMMASSIGN(OUT[3],IN[0]);		\
    }							\
  while(0)

#define ZBW_DER_TMQ_EXP(OUT,IN)				\
  do							\
    {							\
      BI_COLOR_SUMMASSIGN(OUT[0],IN[0]);		\
      BI_COLOR_SUMMASSIGN(OUT[1],IN[1]);		\
      BI_COLOR_ISUBTASSIGN(OUT[2],IN[0]);		\
      BI_COLOR_ISUMMASSIGN(OUT[3],IN[1]);		\
    }							\
  while(0)

