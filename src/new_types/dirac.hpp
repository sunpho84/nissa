#ifndef _DIRAC_HPP
#define _DIRAC_HPP

#include <stdio.h>
#include <string.h>

#include "complex.hpp"
#include "base/debug.hpp"

#ifndef EXTERN_DIRAC
 #define EXTERN_DIRAC extern
 #define ONLY_INSTANTIATION
#endif

#define NDIRAC 4

namespace nissa
{
  //The structure for gamma matrix
  struct dirac_matr
  {
    int pos[NDIRAC];
    complex entr[NDIRAC];
  };
  
  //The base of the 16 gamma matrixes, the two rotators and Ci=G0*Gi*G5
  CUDA_MANAGED EXTERN_DIRAC dirac_matr base_gamma[19];
  EXTERN_DIRAC dirac_matr Pplus,Pminus;
  EXTERN_DIRAC char gtag[19][3]
#ifndef ONLY_INSTANTIATION
  ={"S0","V1","V2","V3","V0","P5","A1","A2","A3","A0","T1","T2","T3","B1","B2","B3","C1","C2","C3"}
#endif
    ;
  CUDA_MANAGED EXTERN_DIRAC int tau3[2]
#ifndef ONLY_INSTANTIATION
  ={-1,+1}
#endif
    ;
  
  //Initialize a dirac matrix with outside entries
  inline dirac_matr init_dirac(int pos0,double rea0,double ima0,int pos1,double rea1,double ima1,int pos2,double rea2,double ima2,int pos3,double rea3,double ima3)
  {
    dirac_matr out;
    
    out.pos[0]=pos0;
    out.pos[1]=pos1;
    out.pos[2]=pos2;
    out.pos[3]=pos3;
    
    out.entr[0][0]=rea0;
    out.entr[1][0]=rea1;
    out.entr[2][0]=rea2;
    out.entr[3][0]=rea3;
    
    out.entr[0][1]=ima0;
    out.entr[1][1]=ima1;
    out.entr[2][1]=ima2;
    out.entr[3][1]=ima3;
    
    return out;
  }
  
  //If the two dirac matrix in1 and in2 have the same position structure sum
  //them, otherwise it crashes
  inline dirac_matr dirac_summ(const dirac_matr& in1,const dirac_matr& in2)
  {
    dirac_matr out;
    
    for(int ig=0;ig<NDIRAC;ig++)
      if(in1.pos[ig]==in2.pos[ig])
	{
	  out.pos[ig]=in1.pos[ig];
	  complex_summ(out.entr[ig],in1.entr[ig],in2.entr[ig]);
	}
      else 
	crash("The two matrix passed to sum have different positions");
    
    return out;
  }
  
  inline dirac_matr dirac_subt(const dirac_matr& in1,const dirac_matr& in2)
  {
    dirac_matr out;
    
    for(int ig=0;ig<NDIRAC;ig++)
      if(in1.pos[ig]==in2.pos[ig])
	{
	  out.pos[ig]=in1.pos[ig];
	  complex_subt(out.entr[ig],in1.entr[ig],in2.entr[ig]);
	}
      else
	crash("The two matrix passed to sum have different positions");
    
    return out;
  }
  
  //Assign to the first dirac matrixes the product of the second and the third
  CUDA_HOST_AND_DEVICE inline dirac_matr dirac_prod(const dirac_matr& in1,const dirac_matr& in2)
  {
    dirac_matr out;
    
    //This is the line on the first matrix
    for(int ig1=0;ig1<NDIRAC;ig1++)
      {
	//This is the line to be taken on the second matrix
	int ig2=in1.pos[ig1];
	
	//For each line, the column of the output matrix which is
	//different from 0 is the column of the second matrix different
	//from 0 on the line with index equal to the column of the first
	//matrix which is different from 0 (that is, ig2)
	out.pos[ig1]=in2.pos[ig2];
	
	//The entries of the output is, on each line, the complex
	//product of the entries of the first matrix on that line, for
	//the entries of the second matrix on the line with the index
	//equal to the column of the first matrix which is different
	//from 0 (which again is ig2)
	unsafe_complex_prod(out.entr[ig1],in1.entr[ig1],in2.entr[ig2]);
      }
    
    return out;
  }
  
  inline dirac_matr operator*(const dirac_matr& in1,const dirac_matr& in2)
  {
    return dirac_prod(in1,in2);
  }
  
  inline dirac_matr dirac_prod_double(const dirac_matr& in1,const double& in2)
  {
    dirac_matr out;
    
    for(int id=0;id<NDIRAC;id++)
      {
	out.pos[id]=in1.pos[id];
	complex_prod_double(out.entr[id],in1.entr[id],in2);
      }
    
    return out;
  }
  
  inline dirac_matr dirac_prod_idouble(const dirac_matr& in1,const double& in2)
  {
    dirac_matr out;
    
    for(int id=0;id<NDIRAC;id++)
      {
	out.pos[id]=in1.pos[id];
	complex_prod_idouble(out.entr[id],in1.entr[id],in2);
      }
    
    return out;
  }
  
  inline dirac_matr dirac_prod_complex(const dirac_matr& in1,const complex& in2)
  {
    dirac_matr out;
    
    for(int id=0;id<NDIRAC;id++)
      {
	out.pos[id]=in1.pos[id];
	unsafe_complex_prod(out.entr[id],in1.entr[id],in2);
      }
    
    return out;
  }
  
  //take the hermitian
  inline dirac_matr dirac_herm(const dirac_matr& in)
  {
    dirac_matr out;
    
    for(int id=0;id<NDIRAC;id++)
      {
	for(int jd=id+1;jd<NDIRAC;jd++)
	  if(in.pos[id]==in.pos[jd])
	    crash("pos[%d]=%d==pos[%d]",id,in.pos[id],jd);
	int od=in.pos[id];
	out.pos[od]=id;
	complex_conj(out.entr[od],in.entr[id]);
      }
    
    return out;
  }
  
  //Assign to the first dirac the product of the second by the complex
  //number passed as argument
  inline dirac_matr dirac_compl_prod(const dirac_matr& in,const complex& c)
  {
    dirac_matr out;
    
    //This is the line on the matrix
    for(int ig=0;ig<NDIRAC;ig++)
      {
	out.pos[ig]=in.pos[ig];
	
	safe_complex_prod(out.entr[ig],in.entr[ig],c);
      }
    
    return out;
  }
  
  //Print the dirac matrix passed as argument only on node 0
  inline void print_dirac(const dirac_matr& in)
  {
    for(int ir=0;ir<NDIRAC;ir++)
      {
	int pos=in.pos[ir];
	
	for(int ic=0;ic<pos;ic++)
	  printf("+%02.2f,+%02.2f\t",0.,0.);
	
	printf("%+02.2f,%+02.2f\t",in.entr[ir][0],in.entr[ir][1]);
	
	for(int ic=pos+1;ic<NDIRAC;ic++)
	  printf("+%02.2f,+%02.2f\t",0.,0.);
	
	printf("\n");
      }
  }
  
  void init_base_gamma();
  
  /// Pair of dirac matrices
  struct idirac_pair_t
  {
    int si;
    int so;
    idirac_pair_t(int si,int so) : si(si),so(so) {}
  };
}

#undef EXTERN_DIRAC

#endif
