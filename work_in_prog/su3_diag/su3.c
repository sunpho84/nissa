#include "nissa.h"

void color_print(color i)
{
  for(int j=0;j<3;j++)
    printf("%d: %lg %lg\n",j,i[j][0],i[j][1]);
}

/*
void su3_maximal_trace_project_iteration(su3 o,su3 i)
{
  //compute the first two entries of A0
  complex  a={ i[0][0][0]+i[1][1][0],-i[0][0][1]+i[1][1][1]};
  complex  b={-i[0][1][0]+i[1][0][0],-i[0][1][1]-i[1][0][1]};
  double IN0=1/sqrt(a[0]*a[0]+a[1]*a[1]+b[0]*b[0]+b[1]*b[1]);
  a[0]*=IN0;a[1]*=IN0;b[0]*=IN0;b[1]*=IN0;
  //compute the first two entries of A1
  complex c ={ i[1][1][0]+i[2][2][0],-i[1][1][1]+i[2][2][1]};
  complex d ={-i[1][2][0]+i[2][1][0],-i[1][2][1]-i[2][1][1]};
  double IN1=1/sqrt(c[0]*c[0]+c[1]*c[1]+d[0]*d[0]+d[1]*d[1]);
  c[0]*=IN1;c[1]*=IN1;d[0]*=IN1;d[1]*=IN1;
  //compute the first two entries of A2
  complex  e={ i[0][0][0]+i[2][2][0],-i[0][0][1]+i[2][2][1]};
  complex  f={-i[0][2][0]+i[2][0][0],-i[0][2][1]-i[2][0][1]};
  double IN2=1/sqrt(e[0]*e[0]+e[1]*e[1]+f[0]*f[0]+f[1]*f[1]);
  e[0]*=IN2;e[1]*=IN2;f[0]*=IN2;f[1]*=IN2;
  //compute df* and de* which are needed subsequentially
  complex dfc,dec;
  unsafe_complex_conj2_prod(dfc,d,f);
  unsafe_complex_conj2_prod(dec,d,e);
  //now just fill the output matrix
  unsafe_complex_prod(o[0][0],a,e);complex_subt_the_prod(o[0][0],b,dfc);
  unsafe_complex_prod(o[0][1],b,c);
  unsafe_complex_prod(o[0][2],a,f);complex_summ_the_prod(o[0][2],b,dec);
  unsafe_complex_conj2_prod_minus(o[1][0],e,b);complex_subt_the_conj2_prod(o[1][0],dfc,a);
  unsafe_complex_conj2_prod(o[1][1],c,a);
  unsafe_complex_conj2_prod_minus(o[1][2],f,b);complex_summ_the_conj2_prod(o[1][2],dec,a);
  unsafe_complex_conj_conj_prod_minus(o[2][0],c,f);
  o[2][1][0]=-d[0];o[2][1][1]=d[1];
  unsafe_complex_conj_conj_prod(o[2][2],c,e);
  
  su3_print(o);
  su3 io;
  safe_su3_hermitian(io,o);
  su3 id;
  su3_prod_su3(id,io,o);
  
  su3_print(id);
  crash("");
  
  //Exponentiate. Commented lines are original from APE (more or less)
  up0=(a00+a11~)
  up1=(a01-a10~)
  icsi=1/sqrt(|a00+a11~|^2+|a01-a10~|^2)
  up0=up0*icsi
  appuno=up0~
  up1=up1*icsi
  appdue=-up1
  n0=a00*appuno+a10*appdue
  n2=a02*appuno+a12*appdue
  n3=-a00*appdue~+a10*appuno~
  n6=a20
  n5=-a02*appdue~+a12*appuno~
  n7=a21
  up1=n5-n7~
  n4=-a01*appdue~+a11*appuno~
  n8=a22
  up0=n8~+n4 
  icsi=1/sqrt(up0*up0~+up1*up1~)
  up0=up0*icsi
  bppuno=up0~
  up1=up1*icsi
  bppdue=-up1
  a02=n2
  a20=-n3*bppdue~+n6*bppuno~
  up1=a02-a20~
  a00=n0
  a22=-n5*bppdue~+n8*bppuno~
  up0=a22~+a00
  icsi=1/sqrt(up0*up0~+up1*up1~)
  up0=up0*icsi
  cppuno=up0~
  up1=up1*icsi
  cppdue=-up1
  e0=appuno
  e6=bppdue~*appdue~
  g00=cppuno*e0+cppdue*e6
}
*/

void diagonalize_su3(color eival,su3 eivec,su3 U)
{
  double ts3=3*sqrt(3);
  complex f={-0.5,-0.5*sqrt(3)};
  
  //find eigenvalues
  complex d={U[0][0][0]+U[1][1][0]+U[2][2][0],U[0][0][1]+U[1][1][1]+U[2][2][1]};
  double dm=d[0]*d[0]+d[1]*d[1];
  complex d2;unsafe_complex_prod(d2,d,d);
  complex d3;unsafe_complex_prod(d3,d2,d);
  
  complex arg={27-dm*(18+dm)+8*d3[0],0};
  complex_sqrt(arg,arg);
  complex a={(-27+9*dm-2*d3[0]+ts3*arg[0])/2,(-2*d3[1]+ts3*arg[1])/2};
  complex_pow(a,a,1.0/3);
  complex b;
  complex_reciprocal(b,a);
  complex c={-d2[0]+3*d[0],-d2[1]-3*d[1]};
  safe_complex_prod(b,b,c);
  
  //this is the first solution.
  //to get other three combine d b a differently
  eival[0][0]=(d[0]+b[0]-a[0])/3;
  eival[0][1]=(d[1]+b[1]-a[1])/3;
  
  safe_complex_prod(a,a,f);
  safe_complex_conj2_prod(b,b,f);
  eival[1][0]=(d[0]+b[0]-a[0])/3;
  eival[1][1]=(d[1]+b[1]-a[1])/3;
  
  safe_complex_prod(a,a,f);
  safe_complex_conj2_prod(b,b,f);
  eival[2][0]=(d[0]+b[0]-a[0])/3;
  eival[2][1]=(d[1]+b[1]-a[1])/3;
  
  //check
  /*
  for(int i=0;i<3;i++)
    {
      complex res={1,0};
      complex_subt_the_conj1_prod(res,d,eival[i]);
      complex temp;
      unsafe_complex_prod(temp,eival[i],eival[i]);
      complex_summ_the_prod(res,d,temp);
      complex_subt_the_prod(res,temp,eival[i]);
      printf("Res: %lg %lg\n",res[0],res[1]);
    }
  */
  if(eivec!=NULL)
    for(int i=0;i<3;i++)
      {
	eivec[i][i][0]=-1;
	eivec[i][i][1]= 0;
	
	int j=(i+1)%3;
	int k=(j+1)%3;
	
	complex Ujj,Ukk;
	complex_subt(Ujj,U[j][j],eival[i]);
	complex_subt(Ukk,U[k][k],eival[i]);
	
	//compute the three determinants
	complex D;
	unsafe_complex_prod(  D,Ujj,Ukk);
	complex_subt_the_prod(D,U[j][k],U[k][j]);
	complex Dj;
	unsafe_complex_prod(  Dj,U[j][i],Ukk);
	complex_subt_the_prod(Dj,U[j][k],U[k][i]);
	complex Dk;
	unsafe_complex_prod(  Dk,Ujj,U[k][i]);
	complex_subt_the_prod(Dk,U[j][i],U[k][j]);
	
	//find the roots
	complex_reciprocal( eivec[i][k],D);
	unsafe_complex_prod(eivec[i][j],eivec[i][k],Dj);
	safe_complex_prod(  eivec[i][k],eivec[i][k],Dk);
	
	//normalize the eigenvector to 1
	double inv_norm=1/sqrt(1+eivec[i][j][0]*eivec[i][j][0]+eivec[i][j][1]*eivec[i][j][1]
			       + eivec[i][k][0]*eivec[i][k][0]+eivec[i][k][1]*eivec[i][k][1]);
	complex_prod_real(eivec[i][0],eivec[i][0],inv_norm);
	complex_prod_real(eivec[i][1],eivec[i][1],inv_norm);
	complex_prod_real(eivec[i][2],eivec[i][2],inv_norm);
	
	//check that this is truly the eigenvector
	/*
	color test1,test2,diff;
	unsafe_su3_prod_color(test1,U,eivec[i]);
	unsafe_color_prod_complex(test2,eivec[i],eival[i]);
	color_subt(diff,test1,test2);
	printf("Testing %d, eig: %lg %lg\n",i,eival[i][0],eival[i][1]);
	//color_print(eivec[i]);
	color_print(diff);
	*/
      }
}

void reconstruct_su3(su3 U,color eival,su3 eivec)
{
  su3 temp;
  for(int i=0;i<3;i++)
    unsafe_color_prod_complex(temp[i],eivec[i],eival[i]);
  su3_dag_prod_su3(U,eivec,temp);
}

int main(int narg,char **arg)
{
  init_nissa();
  
  //set the lattice grid
  int T=8;
  int L=4;
  
  //init the grid
  init_grid(T,L);

  quad_su3 *conf=nissa_malloc("conf",loc_vol+bord_vol+edge_vol,quad_su3);
  
  //read conf, compute plaquette, print it
  read_ildg_gauge_conf(conf,"../../data/L4T8conf");
  
  color eival;
  su3 eivec;
  diagonalize_su3(eival,eivec,conf[0][0]);
  
  su3 U;
  reconstruct_su3(U,eival,eivec);
  //su3_subt(U,U,conf[0][0]);
  su3_print(U);
  su3_print(conf[0][0]);
  
  nissa_free(conf);
  
  close_nissa();
  
  return 0;
}
