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
  
  /*
  su3_print(o);
  su3 io;
  safe_su3_hermitian(io,o);
  su3 id;
  su3_prod_su3(id,io,o);
  
  su3_print(id);
  crash("");
  */
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

void diagonalize_su3(su3 U)
{
  su3_print(U);
  double ts3=3*sqrt(3);
  
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
  complex sol1={(d[0]+b[0]-a[0])/3,(d[1]+b[1]-a[1])/3};
  printf("%lg %lg\n",sol1[0],sol1[1]);
}

