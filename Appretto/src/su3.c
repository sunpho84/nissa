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

//////////////////////////////////// Put to zero /////////////////////////////////

void color_put_to_zero(color m){memset(m,0,sizeof(color));}
void su3_put_to_zero(su3 m){memset(m,0,sizeof(su3));}
void as2t_su3_put_to_zero(as2t_su3 m){memset(m,0,sizeof(as2t_su3));}
void spincolor_put_to_zero(spincolor m){memset(m,0,sizeof(spincolor));}
void su3spinspin_put_to_zero(su3spinspin m){memset(m,0,sizeof(su3spinspin));}

//////////////////////////////////////// Copy /////////////////////////////////////

void color_copy(color b,color a){memcpy(b,a,sizeof(color));}
void su3_copy(su3 b,su3 a){memcpy(b,a,sizeof(su3));}
void quad_su3_copy(quad_su3 b,quad_su3 a){memcpy(b,a,sizeof(quad_su3));}
void spincolor_copy(spincolor b,spincolor a){memcpy(b,a,sizeof(spincolor));}

////////////////////////////////// Operations between colors //////////////////////////

//summ two colors
void color_summ(color a,color b,color c)
{for(int i=0;i<6;i++) ((double*)a)[i]=((double*)b)[i]+((double*)c)[i];}

void color_isumm(color a,color b,color c)
{for(int i=0;i<6;i+=2) {((double*)a)[i]=((double*)b)[i]-((double*)c)[i+1];((double*)a)[i+1]=((double*)b)[i+1]+((double*)c)[i];}}

void color_isubt(color a,color b,color c)
{for(int i=0;i<6;i+=2) {((double*)a)[i]=((double*)b)[i]+((double*)c)[i+1];((double*)a)[i+1]=((double*)b)[i+1]-((double*)c)[i];}}

void color_subt(color a,color b,color c)
{for(int i=0;i<6;i++) ((double*)a)[i]=((double*)b)[i]-((double*)c)[i];}

void summassign_color(color a,color b)
{for(int i=0;i<6;i++) ((double*)a)[i]+=((double*)b)[i];}

void subtassign_color(color a,color b)
{for(int i=0;i<6;i++) ((double*)a)[i]-=((double*)b)[i];}

void summassign_icolor(color a,color b)
{for(int i=0;i<6;i+=2) {((double*)a)[i]-=((double*)b)[i+1];((double*)a)[i+1]+=((double*)b)[i];}}

void subtassign_icolor(color a,color b)
{for(int i=0;i<6;i+=2) {((double*)a)[i]+=((double*)b)[i+1];((double*)a)[i+1]-=((double*)b)[i];}}

////////////////////////////////// Operations between su3 //////////////////////////

//summ two su3 matrixes
void su3_summ(su3 a,su3 b,su3 c)
{for(int i=0;i<18;i++) ((double*)a)[i]=((double*)b)[i]+((double*)c)[i];}

//Product of two su3 matrixes
void su3_prod_su3(su3 a,su3 b,su3 c)
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
void su3_dag_prod_su3(su3 a,su3 b,su3 c)
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
void su3_prod_su3_dag(su3 a,su3 b,su3 c)
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
void su3_dag_prod_su3_dag(su3 a,su3 b,su3 c)
{
  for(int ir_out=0;ir_out<3;ir_out++)
    for(int ic_out=0;ic_out<3;ic_out++)
      {
	unsafe_complex_conj_conj_prod(a[ir_out][ic_out],b[0][ir_out],c[ic_out][0]);
	for(int itemp=1;itemp<3;itemp++)
	  complex_summ_the_conj_conj_prod(a[ir_out][ic_out],b[itemp][ir_out],c[ic_out][itemp]);
      }
}

//product of an su3 matrix by a complex
void safe_su3_prod_complex(su3 a,su3 b,complex c)
{
  complex *ca=(complex*)a;
  complex *cb=(complex*)b;

  for(int i=0;i<9;i++) safe_complex_prod(ca[i],cb[i],c);
}

//product of an su3 matrix by a real
void su3_prod_real(su3 a,su3 b,double r)
{
  double *da=(double*)a;
  double *db=(double*)b;

  for(int i=0;i<18;i++) da[i]=db[i]*r;
}

////////////////////// products between su3 and color //////////////////

//product of an su3 matrix by a color vector
void unsafe_su3_prod_color(color a,su3 b,color c)
{
  for(int c1=0;c1<3;c1++)
    {
      unsafe_complex_prod(a[c1],b[c1][0],c[0]);
      for(int c2=1;c2<3;c2++) complex_summ_the_prod(a[c1],b[c1][c2],c[c2]);
    }
}

void unsafe_color_prod_su3(color a,color b,su3 c)
{
  for(int c1=0;c1<3;c1++)
    {
      unsafe_complex_prod(a[c1],b[0],c[0][c1]);
      for(int c2=1;c2<3;c2++) complex_summ_the_prod(a[c1],b[c2],c[c2][c1]);
    }
}

void unsafe_summ_su3_prod_color(color a,su3 b,color c)
{for(int c1=0;c1<3;c1++) for(int c2=0;c2<3;c2++) complex_summ_the_prod(a[c1],b[c1][c2],c[c2]);}

void unsafe_subt_su3_prod_color(color a,su3 b,color c)
{for(int c1=0;c1<3;c1++) for(int c2=0;c2<3;c2++) complex_subt_the_prod(a[c1],b[c1][c2],c[c2]);}

//product of an su3 matrix by a color vector
void safe_su3_prod_color(color a,su3 b,color c)
{
  color t;
  for(int c1=0;c1<3;c1++)
    {
      unsafe_complex_prod(t[c1],b[c1][0],c[0]);
      for(int c2=1;c2<3;c2++) complex_summ_the_prod(t[c1],b[c1][c2],c[c2]);
    }
  color_copy(a,t);
}

//product of an su3 matrix by a color vector
void unsafe_su3_dag_prod_color(color a,su3 b,color c)
{
  for(int c1=0;c1<3;c1++)
    {
      unsafe_complex_conj1_prod(a[c1],b[0][c1],c[0]);
      for(int c2=1;c2<3;c2++) complex_summ_the_conj1_prod(a[c1],b[c2][c1],c[c2]);
    }
}

void unsafe_color_prod_su3_dag(color a,color b,su3 c)
{
  for(int c1=0;c1<3;c1++)
    {
      unsafe_complex_conj2_prod(a[c1],b[0],c[c1][0]);
      for(int c2=1;c2<3;c2++) complex_summ_the_conj2_prod(a[c1],b[c2],c[c1][c2]);
    }
}

void unsafe_summ_su3_dag_prod_color(color a,su3 b,color c)
{for(int c1=0;c1<3;c1++) for(int c2=0;c2<3;c2++) complex_summ_the_conj1_prod(a[c1],b[c2][c1],c[c2]);}

//product of an su3 matrix by a color vector
void safe_su3_dag_prod_color(color a,su3 b,color c)
{
  color t;
  for(int c1=0;c1<3;c1++)
    {
      safe_complex_conj1_prod(t[c1],b[0][c1],c[0]);
      for(int c2=1;c2<3;c2++) complex_summ_the_conj1_prod(t[c1],b[c2][c1],c[c2]);
    }
  color_copy(a,t);
}

////////////////////////////////// Operations between spincolor ///////////////////

//summ two spincolors
void spincolor_summ(spincolor a,spincolor b,spincolor c)
{for(int i=0;i<24;i++) ((double*)a)[i]=((double*)b)[i]+((double*)c)[i];}

//subtract two spincolors
void spincolor_subt(spincolor a,spincolor b,spincolor c)
{for(int i=0;i<24;i++) ((double*)a)[i]=((double*)b)[i]-((double*)c)[i];}

//summ after multyplying for a complex factor
void safe_spincolor_summ_with_cfactor(spincolor a,spincolor b,spincolor c,complex factor)
{
  complex temp;
  for(int i=0;i<12;i++)
    {
      unsafe_complex_prod(temp,factor,((complex*)c)[i]);
      complex_summ(((complex*)a)[i],((complex*)b)[i],temp);
    }
}

//summ after multyplying for a real factor
void safe_spincolor_summ_with_rfactor(spincolor a,spincolor b,spincolor c,double factor)
{for(int i=0;i<24;i++) ((double*)a)[i]=factor*((double*)b)[i]+((double*)c)[i];}

//spincolor*complex
void safe_complex_prod_spincolor(spincolor out,complex a,spincolor in)
{for(int i=0;i<12;i++) safe_complex_prod(((complex*)out)[i],a,((complex*)in)[i]);}

//spincolor*real
void assign_spincolor_prod_real(spincolor out,double factor)
{for(int i=0;i<24;i++) ((double*)out)[i]*=factor;}

//summ assign
void summassign_spincolor(spincolor out,spincolor in)
{for(int i=0;i<24;i++) ((double*)out)[i]+=((double*)in)[i];}

//subt assign
void subtassign_spincolor(spincolor out,spincolor in)
{for(int i=0;i<24;i++) ((double*)out)[i]-=((double*)in)[i];}

//spincolor*real
void unsafe_summassign_spincolor_prod_real(spincolor out,spincolor in,double factor)
{for(int i=0;i<24;i++) ((double*)out)[i]+=((double*)in)[i]*factor;}

//spincolor*i*real
void unsafe_summassign_spincolor_prod_ireal(spincolor out,spincolor in,double factor)
{
  for(int i=0;i<24;i+=2)
    {
      ((double*)out)[i  ]-=((double*)in)[i+1]*factor;
      ((double*)out)[i+1]+=((double*)in)[i  ]*factor;
    }
}

//spincolor*i*real
void unsafe_spincolor_summ_with_ifactor(spincolor out,spincolor a,spincolor b,double factor)
{
  for(int i=0;i<24;i+=2)
    {
      ((double*)out)[i  ]=((double*)a)[i  ]-((double*)b)[i+1]*factor;
      ((double*)out)[i+1]=((double*)a)[i+1]+((double*)b)[i  ]*factor;
    }
}

//dirac*spincolor
void unsafe_dirac_prod_spincolor(spincolor out,dirac_matr *m,spincolor in)
{for(int id1=0;id1<4;id1++) for(int ic=0;ic<3;ic++) safe_complex_prod(out[id1][ic],m->entr[id1],in[m->pos[id1]][ic]);}

//dirac*spincolor
void safe_dirac_prod_spincolor(spincolor out,dirac_matr *m,spincolor in)
{
  spincolor tmp;
  unsafe_dirac_prod_spincolor(tmp,m,in);
  spincolor_copy(out,tmp);
}

//su3*spincolor
void unsafe_su3_prod_spincolor(spincolor out,su3 U,spincolor in)
{for(int is=0;is<4;is++) unsafe_su3_prod_color(out[is],U,in[is]);}

void unsafe_summ_su3_prod_spincolor(spincolor out,su3 U,spincolor in)
{for(int is=0;is<4;is++) unsafe_summ_su3_prod_color(out[is],U,in[is]);}

//su3^*spincolor
void unsafe_su3_dag_prod_spincolor(spincolor out,su3 U,spincolor in)
{for(int is=0;is<4;is++) unsafe_su3_dag_prod_color(out[is],U,in[is]);}

void unsafe_summ_su3_dag_prod_spincolor(spincolor out,su3 U,spincolor in)
{for(int is=0;is<4;is++) unsafe_summ_su3_dag_prod_color(out[is],U,in[is]);}

//su3^*gamma*spincolor
void unsafe_su3_dag_dirac_prod_spincolor(spincolor out,su3 U,dirac_matr *m,spincolor in)
{
  color tmp;
  for(int id1=0;id1<4;id1++)
    {
      for(int ic=0;ic<3;ic++) unsafe_complex_prod(tmp[ic],m->entr[id1],in[m->pos[id1]][ic]);
      unsafe_su3_dag_prod_color(out[id1],U,tmp);
    }
}

void unsafe_summ_su3_dag_dirac_prod_spincolor(spincolor out,su3 U,dirac_matr *m,spincolor in)
{
  color tmp;
  for(int id1=0;id1<4;id1++)
    {
      for(int ic=0;ic<3;ic++) unsafe_complex_prod(tmp[ic],m->entr[id1],in[m->pos[id1]][ic]);
      unsafe_summ_su3_dag_prod_color(out[id1],U,tmp);
    }
}

//su3*dirac*spincolor
void unsafe_su3_dirac_prod_spincolor(spincolor out,su3 U,dirac_matr *m,spincolor in)
{
  color tmp;
  for(int id1=0;id1<4;id1++)
    {
      for(int ic=0;ic<3;ic++) unsafe_complex_prod(tmp[ic],m->entr[id1],in[m->pos[id1]][ic]);
      unsafe_su3_prod_color(out[id1],U,tmp);
    }
}

void unsafe_subt_su3_dirac_prod_spincolor(spincolor out,su3 U,dirac_matr *m,spincolor in)
{
  color tmp;
  for(int id1=0;id1<4;id1++)
    {
      for(int ic=0;ic<3;ic++) unsafe_complex_prod(tmp[ic],m->entr[id1],in[m->pos[id1]][ic]);
      unsafe_subt_su3_prod_color(out[id1],U,tmp);
    }
}

/////////////////////////////////////// Complicated things /////////////////////

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

  su3_prod_su3(AB,conf[A][mu],conf[B][nu]);
  su3_prod_su3(AC,conf[A][nu],conf[C][mu]);
  su3_prod_su3_dag(square,AB,AC);
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
	      su3_prod_su3(temp1,conf[X][mu],conf[A][nu]);         //    B--<--Y 
	      su3_prod_su3_dag(temp2,temp1,conf[B][mu]);           //    |  1  | 
	      su3_prod_su3_dag(temp1,temp2,conf[X][nu]);           //    |     | 
	      su3_summ(leaves_summ,leaves_summ,temp1);	           //    X-->--A 

	      //Leaf 2
	      su3_prod_su3_dag(temp1,conf[X][nu],conf[C][mu]);      //   C--<--B
	      su3_prod_su3_dag(temp2,temp1,conf[D][nu]);            //   |  2  | 
	      su3_prod_su3(temp1,temp2,conf[D][mu]);		    //   |     | 
	      su3_summ(leaves_summ,leaves_summ,temp1);		    //   D-->--X
	      
	      //Leaf 3
	      su3_dag_prod_su3_dag(temp1,conf[D][mu],conf[E][nu]);  //   D--<--X
	      su3_prod_su3(temp2,temp1,conf[E][mu]);		    //   |  3  | 
	      su3_prod_su3(temp1,temp2,conf[F][nu]);		    //   |     | 
	      su3_summ(leaves_summ,leaves_summ,temp1);		    //   E-->--F
	      
	      //Leaf 4
	      su3_dag_prod_su3(temp1,conf[F][nu],conf[F][mu]);       //  X--<--A 
	      su3_prod_su3(temp2,temp1,conf[G][nu]);                 //  |  4  | 
	      su3_prod_su3_dag(temp1,temp2,conf[X][mu]);             //  |     |  
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
	  unsafe_su3_prod_color(temp_d1,Pmunu[imunu],in[smunu_pos[d1][imunu]]);
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
