#ifndef _DIRAC_HPP
#define _DIRAC_HPP

#include "complex.hpp"
#include "spin.hpp"
#include "su3.hpp"

namespace nissa
{
  //Initialize a dirac matrix with outside entries
  inline void init_dirac(dirac_matr *out,int pos0,double rea0,double ima0,int pos1,double rea1,double ima1,int pos2,double rea2,double ima2,int pos3,double rea3,double ima3)
  {
    out->pos[0]=pos0;
    out->pos[1]=pos1;
    out->pos[2]=pos2;
    out->pos[3]=pos3;
    
    out->entr[0][0]=rea0;
    out->entr[1][0]=rea1;
    out->entr[2][0]=rea2;
    out->entr[3][0]=rea3;
    
    out->entr[0][1]=ima0;
    out->entr[1][1]=ima1;
    out->entr[2][1]=ima2;
    out->entr[3][1]=ima3;
  }
  
  //If the two dirac matrix in1 and in2 have the same position structure sum
  //them, otherwise it crashes
  inline void dirac_summ(dirac_matr *out,dirac_matr *in1,dirac_matr *in2)
  {
    for(int ig=0;ig<4;ig++)
      if(in1->pos[ig]==in2->pos[ig])
	{
	  out->pos[ig]=in1->pos[ig];
	  complex_summ(out->entr[ig],in1->entr[ig],in2->entr[ig]);
	}
      else 
	crash("The two matrix passed to sum have different positions");
  }
  inline void dirac_subt(dirac_matr *out,dirac_matr *in1,dirac_matr *in2)

    
  {
    for(int ig=0;ig<4;ig++)
      if(in1->pos[ig]==in2->pos[ig])
	{
	  out->pos[ig]=in1->pos[ig];
	  complex_subt(out->entr[ig],in1->entr[ig],in2->entr[ig]);
	}
      else 
	crash("The two matrix passed to sum have different positions");
  }
  
  //Assign to the first dirac matrixes the product of the second and the third
  inline void dirac_prod(dirac_matr *out,dirac_matr *in1,dirac_matr *in2)
  {
    dirac_matr temp; //this is needed to avoid to overwrite one of the input
    
    //This is the line on the first matrix
    for(int ig1=0;ig1<4;ig1++)
      {
	//This is the line to be taken on the second matrix
	int ig2=in1->pos[ig1];
	
	//For each line, the column of the output matrix which is
	//different from 0 is the column of the second matrix different
	//from 0 on the line with index equal to the column of the first
	//matrix which is different from 0 (that is, ig2)
	temp.pos[ig1]=in2->pos[ig2];
	
	//The entries of the output is, on each line, the complex
	//product of the entries of the first matrix on that line, for
	//the entries of the second matrix on the line with the index
	//equal to the column of the first matrix which is different
	//from 0 (which again is ig2)
	unsafe_complex_prod(temp.entr[ig1],in1->entr[ig1],in2->entr[ig2]);
      }
    
    memcpy(out->pos,temp.pos,sizeof(int)*4);
    memcpy(out->entr,temp.entr,sizeof(complex)*4);  
  }
  
  //Assign to the first dirac the product of the second by the complex
  //number passed as argument
  inline void safe_dirac_compl_prod(dirac_matr *out,dirac_matr *in,complex c)
  {
    //This is the line on the matrix
    for(int ig=0;ig<4;ig++)
      {
	out->pos[ig]=in->pos[ig];
	
	safe_complex_prod(out->entr[ig],in->entr[ig],c);
      }
  }
  
  //Assign to the first dirac the product of the second by the complex
  //number passed as argument
  inline void unsafe_dirac_compl_prod(dirac_matr *out,dirac_matr *in,complex c)
  {
    //This is the line on the matrix
    for(int ig=0;ig<4;ig++)
      {
	out->pos[ig]=in->pos[ig];
	
	unsafe_complex_prod(out->entr[ig],in->entr[ig],c);
      }
  }
  
  //Summ the passed gamma multiplied by a double to spinspin
  inline void spinspin_dirac_summ_the_prod_double(spinspin out,dirac_matr *in,double r)
  {
    //This is the line on the matrix
    for(int ig=0;ig<4;ig++)
      complex_summ_the_prod_double(out[ig][in->pos[ig]],in->entr[ig],r);
  }
  inline void spinspin_dirac_summ_the_prod_idouble(spinspin out,dirac_matr *in,double r)
  {
    for(int ig=0;ig<4;ig++)
      complex_summ_the_prod_idouble(out[ig][in->pos[ig]],in->entr[ig],r);
  }
  inline void spinspin_dirac_summ_the_prod_complex(spinspin out,dirac_matr *in,complex c)
  {
    for(int ig=0;ig<4;ig++)
      complex_summ_the_prod(out[ig][in->pos[ig]],in->entr[ig],c);
  }
  inline void spinspin_dirac_subt_the_prod_complex(spinspin out,dirac_matr *in,complex c)
  {
    for(int ig=0;ig<4;ig++)
      complex_subt_the_prod(out[ig][in->pos[ig]],in->entr[ig],c);
  }
  inline void spinspin_dirac_prod_double(spinspin out,dirac_matr *in,double r)
  {spinspin_put_to_zero(out);spinspin_dirac_summ_the_prod_double(out,in,r);}
  inline void spinspin_dirac_prod_idouble(spinspin out,dirac_matr *in,double r)
  {spinspin_put_to_zero(out);spinspin_dirac_summ_the_prod_idouble(out,in,r);}
  inline void spinspin_dirac_prod_complex(spinspin out,dirac_matr *in,complex c)
  {spinspin_put_to_zero(out);spinspin_dirac_summ_the_prod_complex(out,in,c);}
  
  //out=m*in
  inline void unsafe_dirac_prod_spinspin(spinspin out,dirac_matr *m,spinspin in)
  {
    for(int id1=0;id1<4;id1++)
      for(int id2=0;id2<4;id2++)
	unsafe_complex_prod(out[id1][id2],m->entr[id1],in[m->pos[id1]][id2]);
  }
  inline void safe_dirac_prod_spinspin(spinspin out,dirac_matr *m,spinspin in)
  {spinspin temp;unsafe_dirac_prod_spinspin(temp,m,in);spinspin_copy(out,temp);}
  inline void unsafe_dirac_prod_colorspinspin(colorspinspin out,dirac_matr *m,colorspinspin in)
  {for(int ic=0;ic<3;ic++) unsafe_dirac_prod_spinspin(out[ic],m,in[ic]);}
  inline void safe_dirac_prod_colorspinspin(colorspinspin out,dirac_matr *m,colorspinspin in)
  {colorspinspin temp;unsafe_dirac_prod_colorspinspin(temp,m,in);colorspinspin_copy(out,temp);}
  inline void unsafe_dirac_prod_su3spinspin(su3spinspin out,dirac_matr *m,su3spinspin in)
  {for(int ic=0;ic<3;ic++) unsafe_dirac_prod_colorspinspin(out[ic],m,in[ic]);}
  inline void safe_dirac_prod_su3spinspin(su3spinspin out,dirac_matr *m,su3spinspin in)
  {su3spinspin temp;unsafe_dirac_prod_su3spinspin(temp,m,in);su3spinspin_copy(out,temp);}

  //out+=m*in
  inline void spinspin_dirac_summ_the_prod_spinspin(spinspin out,dirac_matr *in,spinspin c)
  {
    spinspin temp;
    unsafe_dirac_prod_spinspin(temp,in,c);
    spinspin_summassign(out,temp);
  }

  
  //out=m*in^t
  inline void unsafe_dirac_prod_spinspin_transp(spinspin out,dirac_matr *m,spinspin in)
  {
    for(int id1=0;id1<4;id1++)
      for(int id2=0;id2<4;id2++)
	unsafe_complex_prod(out[id1][id2],m->entr[id1],in[id2][m->pos[id1]]);
  }
  inline void safe_dirac_prod_spinspin_transp(spinspin out,dirac_matr *m,spinspin in)
  {spinspin temp;unsafe_dirac_prod_spinspin_transp(temp,m,in);spinspin_copy(out,temp);}
  
  //out=m*in^+
  inline void unsafe_dirac_prod_spinspin_dag(spinspin out,dirac_matr *m,spinspin in)
  {
    for(int id1=0;id1<4;id1++)
      for(int id2=0;id2<4;id2++)
	unsafe_complex_conj2_prod(out[id1][id2],m->entr[id1],in[id2][m->pos[id1]]);
  }
  inline void safe_dirac_prod_spinspin_dag(spinspin out,dirac_matr *m,spinspin in)
  {spinspin temp;unsafe_dirac_prod_spinspin_dag(temp,m,in);spinspin_copy(out,temp);}
  
  //out=in*m
  inline void unsafe_spinspin_prod_dirac(spinspin out,spinspin in,dirac_matr *m)
  {
    spinspin_put_to_zero(out);
    for(int id1=0;id1<4;id1++)
      for(int id2=0;id2<4;id2++)
	unsafe_complex_prod(out[id1][m->pos[id2]],in[id1][id2],m->entr[id2]);
  }
  inline void safe_spinspin_prod_dirac(spinspin out,spinspin in,dirac_matr *m)
  {spinspin temp;unsafe_spinspin_prod_dirac(temp,in,m);spinspin_copy(out,temp);}
  
  //Print the dirac marix passed as argument only on node 0
  inline void print_dirac(dirac_matr *in)
  {
    for(int ir=0;ir<4;ir++)
      {
	int pos=in->pos[ir];
	for(int ic=0;ic<pos;ic++) printf("+%02.2f,+%02.2f\t",0.,0.);
	printf("%+02.2f,%+02.2f\t",in->entr[ir][0],in->entr[ir][1]);
	for(int ic=pos+1;ic<4;ic++) printf("+%02.2f,+%02.2f\t",0.,0.);
	printf("\n");
      }
  }
  
  //Trace of the product of two spinspins
  inline void summ_the_trace_prod_spinspins(complex c,spinspin a,spinspin b)
  {
    for(int id1=0;id1<4;id1++)
      for(int id2=0;id2<4;id2++)
	complex_summ_the_prod(c,a[id1][id2],b[id2][id1]);
  }
  inline void trace_prod_spinspins(complex c,spinspin a,spinspin b)
  {
    c[0]=c[1]=0;
    summ_the_trace_prod_spinspins(c,a,b);
  }
  inline void trace_spinspin_with_dirac(complex out,spinspin s,dirac_matr *m)
  {
    complex_put_to_zero(out);
    for(int id1=0;id1<4;id1++)
      complex_summ_the_prod(out,m->entr[id1],s[m->pos[id1]][id1]);
  }
  
  void init_base_gamma();
}

#endif
