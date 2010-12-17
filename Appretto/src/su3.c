#pragma once

#include "global.c"

//return the trace of an su3 matrix
void su3_trace(complex tr,su3 m)
{
  tr[0]=m[0][0][0];
  tr[1]=m[0][0][1];

  complex_summ(tr,tr,m[1][1]);
  complex_summ(tr,tr,m[2][2]);
}

//put to zero an su3 matrix
void su3_put_to_zero(su3 m)
{
  memset(m,0,sizeof(su3));
}

//put to zero a color vector
void color_put_to_zero(color m)
{
  memset(m,0,sizeof(color));
}

//put to zero an su3 anti-simmetric tensor
void as2t_su3_put_to_zero(as2t_su3 m)
{
  memset(m,0,sizeof(as2t_su3));
}

//copy a to b
void su3_copy(su3 b,su3 a)
{
  memcpy(b,a,sizeof(su3));
}

//copy a to b
void quad_su3_copy(quad_su3 b,quad_su3 a)
{
  memcpy(b,a,sizeof(quad_su3));
}

//put to zero an su3 matrix of spinspin
void su3spinspin_put_to_zero(su3spinspin m)
{
  memset(m,0,sizeof(su3spinspin));
}

//summ two su3 matrixes
void su3_summ(su3 a,su3 b,su3 c)
{
  for(int i=0;i<18;i++)
    ((double*)a)[i]=((double*)b)[i]+((double*)c)[i];
}

//Product of two su3 matrixes
void su3_su3_prod(su3 a,su3 b,su3 c)
{
  for(int ir_out=0;ir_out<3;ir_out++)
    for(int ic_out=0;ic_out<3;ic_out++)
      {
	unsafe_complex_prod(a[ir_out][ic_out],b[ir_out][0],c[0][ic_out]);
	for(int itemp=1;itemp<3;itemp++)
	  complex_summ_the_prod(a[ir_out][ic_out],b[ir_out][itemp],c[itemp][ic_out]);
      }
}

//Product of two su3 matrixes
void su3_dag_su3_prod(su3 a,su3 b,su3 c)
{
  for(int ir_out=0;ir_out<3;ir_out++)
    for(int ic_out=0;ic_out<3;ic_out++)
      {
	unsafe_complex_conj1_prod(a[ir_out][ic_out],b[0][ir_out],c[0][ic_out]);
	for(int itemp=1;itemp<3;itemp++)
	  complex_summ_the_conj1_prod(a[ir_out][ic_out],b[itemp][ir_out],c[itemp][ic_out]);
      }
}

//Product of two su3 matrixes
void su3_su3_dag_prod(su3 a,su3 b,su3 c)
{
  for(int ir_out=0;ir_out<3;ir_out++)
    for(int ic_out=0;ic_out<3;ic_out++)
      {
	unsafe_complex_conj2_prod(a[ir_out][ic_out],b[ir_out][0],c[ic_out][0]);
	for(int itemp=1;itemp<3;itemp++)
	  complex_summ_the_conj2_prod(a[ir_out][ic_out],b[ir_out][itemp],c[ic_out][itemp]);
      }
}

//Product of two su3 matrixes
void su3_dag_su3_dag_prod(su3 a,su3 b,su3 c)
{
  for(int ir_out=0;ir_out<3;ir_out++)
    for(int ic_out=0;ic_out<3;ic_out++)
      {
	unsafe_complex_conj_conj_prod(a[ir_out][ic_out],b[0][ir_out],c[ic_out][0]);
	for(int itemp=1;itemp<3;itemp++)
	  complex_summ_the_conj_conj_prod(a[ir_out][ic_out],b[itemp][ir_out],c[ic_out][itemp]);
      }
}

//product of an su3 matrix by a color vector
void unsafe_su3_color_prod(color a,su3 b,color c)
{
  for(int c1=0;c1<3;c1++)
    {
      unsafe_complex_prod(a[c1],b[c1][0],c[0]);
      for(int c2=1;c2<3;c2++) complex_summ_the_prod(a[c1],b[c1][c2],c[c2]);
    }
}

//product of an su3 matrix by a color vector
void unsafe_su3_dag_color_prod(color a,su3 b,color c)
{
  for(int c1=0;c1<3;c1++)
    {
      unsafe_complex_conj1_prod(a[c1],b[0][c1],c[0]);
      for(int c2=1;c2<3;c2++) complex_summ_the_conj1_prod(a[c1],b[c2][c1],c[c2]);
    }
}

//product of an su3 matrix by a complex
void safe_su3_complex_prod(su3 a,su3 b,complex c)
{
  complex *ca=(complex*)a;
  complex *cb=(complex*)b;

  for(int i=0;i<9;i++) safe_complex_prod(ca[i],cb[i],c);
}

//square (the proto-plaquette)
/*
     
  C------D
n |      | 
u |      | 
  A--mu--B
  
The square path P_{mu,nu} is defined as U(A,mu)U(B,nu)U^(C,mu)U^(A,nu)=
=U(A,mu)U(B,nu)(U(A,nu)U^(C,mu))^=U(AB,munu)*U^(AC,numu)
*/

void squared_path(su3 square,quad_su3 *conf,int A,int mu,int nu)
{
  int B=loclx_neighup[A][mu];
  int C=loclx_neighup[A][nu];

  su3 AB,AC;

  su3_su3_prod(AB,conf[A][mu],conf[B][nu]);
  su3_su3_prod(AC,conf[A][nu],conf[C][mu]);
  su3_su3_dag_prod(square,AB,AC);
}

//This calculate the global plaquette. It's not done in a very
//efficient way, but it's ok for our scope.
double global_plaquette(quad_su3 *conf)
{
  su3 square;
  complex pl;
  double totlocplaq=0;
  for(int ivol=0;ivol<loc_vol;ivol++)
    for(int idir=0;idir<4;idir++)
      for(int jdir=idir+1;jdir<4;jdir++)
	{
	  squared_path(square,conf,ivol,idir,jdir);
	  su3_trace(pl,square);
	  totlocplaq+=pl[0]/3;
	}
  
  double totplaq;
  MPI_Reduce(&totlocplaq,&totplaq,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
  
  return totplaq/glb_vol/6;
}
//This will calculate 2*a^2*ig*P_{mu,nu}
/*
  ^                   C--<-- B --<--Y 
  |                   |  2  | |  1  | 
  n                   |     | |     | 
  u                   D-->--\X/-->--A 
  |                   D--<--/X\--<--A 
  -----mu---->        |  3  | |  4  | 
  		      |     | |     | 
		      E-->-- F -->--G 
*/
void Pmunu_term(as2t_su3 *Pmunu,quad_su3 *conf)
{
  int A,B,C,D,E,F,G;
  int munu;

  su3 temp1,temp2,leaves_summ;

  for(int X=0;X<loc_vol;X++)
    {
      as2t_su3_put_to_zero(Pmunu[X]);

      munu=0;
      for(int mu=0;mu<4;mu++)
	{
	  A=loclx_neighup[X][mu];
	  D=loclx_neighdw[X][mu];
	  
	  for(int nu=mu+1;nu<4;nu++)
	    {
	      B=loclx_neighup[X][nu];
	      F=loclx_neighdw[X][nu];
	      
	      C=loclx_neighup[D][nu];
	      E=loclx_neighdw[D][nu];

	      G=loclx_neighdw[A][nu];

	      //Put to 0 the summ of the leaves
	      su3_put_to_zero(leaves_summ);
	      
	      //Leaf 1
	      su3_su3_prod(temp1,conf[X][mu],conf[A][nu]);         //    B--<--Y 
	      su3_su3_dag_prod(temp2,temp1,conf[B][mu]);           //    |  1  | 
	      su3_su3_dag_prod(temp1,temp2,conf[X][nu]);           //    |     | 
	      su3_summ(leaves_summ,leaves_summ,temp1);	           //    X-->--A 

	      //Leaf 2
	      su3_su3_dag_prod(temp1,conf[X][nu],conf[C][mu]);      //   C--<--B
	      su3_su3_dag_prod(temp2,temp1,conf[D][nu]);            //   |  2  | 
	      su3_su3_prod(temp1,temp2,conf[D][mu]);		    //   |     | 
	      su3_summ(leaves_summ,leaves_summ,temp1);		    //   D-->--X
	      
	      //Leaf 3
	      su3_dag_su3_dag_prod(temp1,conf[D][mu],conf[E][nu]);  //   D--<--X
	      su3_su3_prod(temp2,temp1,conf[E][mu]);		    //   |  3  | 
	      su3_su3_prod(temp1,temp2,conf[F][nu]);		    //   |     | 
	      su3_summ(leaves_summ,leaves_summ,temp1);		    //   E-->--F
	      
	      //Leaf 4
	      su3_dag_su3_prod(temp1,conf[F][nu],conf[F][mu]);       //  X--<--A 
	      su3_su3_prod(temp2,temp1,conf[G][nu]);                 //  |  4  | 
	      su3_su3_dag_prod(temp1,temp2,conf[X][mu]);             //  |     |  
	      su3_summ(leaves_summ,leaves_summ,temp1);               //  F-->--G 

	      //calculate U-U^dagger
	      for(int ic1=0;ic1<3;ic1++)
		for(int ic2=0;ic2<3;ic2++)
		  {
		    Pmunu[X][munu][ic1][ic2][0]=(leaves_summ[ic1][ic2][0]-leaves_summ[ic2][ic1][0])/4;
		    Pmunu[X][munu][ic1][ic2][1]=(leaves_summ[ic1][ic2][1]+leaves_summ[ic2][ic1][1])/4;
		  }

	      munu++;
	    }
	}
    }
}

//apply the chromo operator to the passed spinor site by site (not yet fully optimized)
void unsafe_apply_point_chromo_operator_to_spincolor(spincolor out,as2t_su3 Pmunu,spincolor in)
{
  color temp_d1;
  
  for(int d1=0;d1<4;d1++)
    {
      color_put_to_zero(out[d1]);
      for(int imunu=0;imunu<6;imunu++)
	{
	  unsafe_su3_color_prod(temp_d1,Pmunu[imunu],in[smunu_pos[d1][imunu]]);
	  for(int c=0;c<3;c++) complex_summ_the_prod(out[d1][c],smunu_entr[d1][imunu],temp_d1[c]);
	}
    }
}

//apply the chromo operator to the passed spinor to the whole volume
void unsafe_apply_chromo_operator_to_spincolor(spincolor *out,as2t_su3 *Pmunu,spincolor *in)
{
  for(int ivol=0;ivol<loc_vol;ivol++) unsafe_apply_point_chromo_operator_to_spincolor(out[ivol],Pmunu[ivol],in[ivol]);
}

//apply the chromo operator to the passed colorspinspin
//normalization as in ape next
void unsafe_apply_chromo_operator_to_colorspinspin(colorspinspin *out,as2t_su3 *Pmunu,colorspinspin *in)
{
  spincolor temp1,temp2;
  
  for(int loc_site=0;loc_site<loc_vol;loc_site++)
    {
      //Loop over the four source dirac indexes
      for(int id_source=0;id_source<4;id_source++) //dirac index of source
	{
	  //Switch the color_spinspin into the spincolor.
	  get_spincolor_from_colorspinspin(temp1,in[loc_site],id_source);
	  
	  unsafe_apply_point_chromo_operator_to_spincolor(temp2,Pmunu[loc_site],temp1);
	  
	  //Switch back the spincolor into the colorspinspin
	  put_spincolor_into_colorspinspin(out[loc_site],temp2,id_source);
	}
    }
}
