#pragma once

void double_to_float_vec(float *out,double *in,int nel)
{for(int iel=0;iel<nel;iel++) out[iel]=(float)(in[iel]);}

void float_to_double_vec(double *out,float *in,int nel)
{for(int iel=0;iel<nel;iel++) out[iel]=(double)(in[iel]);}

//compress an su3 matrix to 8 parameters
void su3_to_su3c(su3c out,su3d in)
{
  out[0]=(float)(in[0][1][0]); //a2
  out[1]=(float)(in[0][1][1]);
  out[2]=(float)(in[0][2][0]); //a3
  out[3]=(float)(in[0][2][1]);
  out[4]=(float)(atan2(in[0][0][1],in[0][0][0])); //theta_a1
  out[5]=(float)(atan2(in[2][0][1],in[2][0][0])); //theta_a2
  out[6]=(float)(in[1][0][0]); //b1
  out[7]=(float)(in[1][0][1]);
}

//decompress 8 parameters to an su3 matrix
void su3c_to_su3f(su3f out,su3c in)
{
  out[0][1][0]=in[0]; //a2
  out[0][1][1]=in[1];
  out[0][2][0]=in[2]; //a3
  out[0][2][1]=in[3];  
  
  float help=in[0]*in[0]+in[1]*in[1]+in[2]*in[2]+in[3]*in[3];
  float N=sqrt(help);
  float rep_N=1/N;
  
  help=sqrt(1-help); //a1
  out[0][0][0]=help*cos(in[4]);
  out[0][0][1]=help*sin(in[4]);
  
  out[1][0][0]=in[6]; //b1
  out[1][0][1]=in[7];
  
  complef p2={rep_N*out[1][0][0],rep_N*out[1][0][1]}; //p2=b1/N
  
  help=sqrt(1-out[0][0][0]*out[0][0][0]-out[0][0][1]*out[0][0][1]-out[1][0][0]*out[1][0][0]-out[1][0][1]*out[1][0][1]);
  out[2][0][0]=help*cos(in[5]); //c1
  out[2][0][1]=help*sin(in[5]);
  
  complef p1={rep_N*out[2][0][0],-rep_N*out[2][0][1]}; //p1=conj(c1)/N
  
  unsafe_complef_conj2_prod(out[1][1],out[0][1],out[0][0]); //b2
  safe_complef_prod(out[1][1],p2,out[1][1]);
  complef_summ_the_conj2_prod(out[1][1],p1,out[0][2]);
  complef_prod_real(out[1][1],out[1][1],-rep_N);

  unsafe_complef_conj2_prod(out[1][2],out[0][2],out[0][0]); //b3
  safe_complef_prod(out[1][2],p2,out[1][2]);
  complef_subt_the_conj2_prod(out[1][2],p1,out[0][1]);
  complef_prod_real(out[1][2],out[1][2],-rep_N);

  unsafe_complef_conj2_prod(out[2][1],out[0][1],out[0][0]); //c2
  safe_complef_conj1_prod(out[2][1],p1,out[2][1]);
  complef_subt_the_conj_conj_prod(out[2][1],p2,out[0][2]);
  complef_prod_real(out[2][1],out[2][1],-rep_N);

  unsafe_complef_conj2_prod(out[2][2],out[0][2],out[0][0]); //c3
  safe_complef_conj1_prod(out[2][2],p1,out[2][2]);
  complef_summ_the_conj_conj_prod(out[2][2],p2,out[0][1]);
  complef_prod_real(out[2][2],out[2][2],-rep_N);
}
