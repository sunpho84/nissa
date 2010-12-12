#pragma once

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
  memset(m,0,18);
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

//This will calculate the clover term
/*
  ^                   C--<-- B --<--X 
  |                   |  2  | |  1  | 
  n                   |     | |     | 
  u                   D-->--\0/-->--A 
  |                   D--<--/0\--<--A 
  -----mu---->        |  3  | |  4  | 
  		      |     | |     | 
		      E-->-- F -->--G 
*/

void clover_term(su3spinspin *clov,quad_su3 *conf)
{
  int A,B,C,D,E,F,G;

  su3 temp1,temp2,leaves_summ;
	      //Put to 0 the summ of the leaves
	      su3_put_to_zero(leaves_summ);


  for(int O=0;O<loc_vol;O++)
    {
      for(int mu=0;mu<4;mu++)
	{
	  A=loclx_neighup[O][mu];
	  D=loclx_neighdw[O][mu];
	  
	  for(int nu=mu+1;nu<4;nu++)
	    {
	      B=loclx_neighup[O][nu];
	      F=loclx_neighdw[O][nu];
	      
	      C=loclx_neighup[D][nu];
	      E=loclx_neighdw[D][nu];

	      G=loclx_neighdw[A][nu];
	      /*
	      //Leaf 1
	      su3_su3_prod(temp1,conf[O][mu],conf[A][nu]);
	      su3_su3_dag_prod(temp2,temp1,conf[B][mu]);
	      su3_su3_dag_prod(temp1,temp2,conf[O][nu]);
	      su3_summ(leaves_summ,leaves_summ,temp1);
	      

	      //Leaf 2
	      su3_su3_dag_prod(temp1,conf[O][nu],conf[C][mu]);
	      su3_su3_dag_prod(temp2,temp1,conf[D][nu]);
	      su3_su3_prod(temp1,temp2,conf[D][mu]);
	      su3_summ(leaves_summ,leaves_summ,temp1);
	      
	      //Leaf 3
	      su3_dag_su3_dag_prod(temp1,conf[D][mu],conf[E][nu]);
	      su3_su3_prod(temp2,temp1,conf[E][mu]);
	      su3_su3_prod(temp1,temp2,conf[F][nu]);
	      su3_summ(leaves_summ,leaves_summ,temp1);
	      */
	      //Leaf 4
	      su3_dag_su3_prod(temp1,conf[F][nu],conf[F][mu]);
	      su3_su3_prod(temp2,temp1,conf[G][nu]);
	      su3_su3_dag_prod(temp1,temp2,conf[O][mu]);
	      su3_summ(leaves_summ,leaves_summ,temp1);
	      
	    }
	}
    }
  
  double total_leaves_summ,partial_summ=leaves_summ[0][0][0];
  MPI_Reduce(&partial_summ,&total_leaves_summ,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);  
  
  if(rank==0) printf("%.12g\n",total_leaves_summ);
}
