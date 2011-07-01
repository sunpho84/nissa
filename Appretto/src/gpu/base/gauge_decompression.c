//decompress 8 parameters to an su3 matrix
{
  real tmp;
  
  out[0][1][0]=in[0]; //a2
  out[0][1][1]=in[1];
  out[0][2][0]=in[2]; //a3
  out[0][2][1]=in[3];  
  
  real help=in[0]*in[0]+in[1]*in[1]+in[2]*in[2]+in[3]*in[3];
  real rep_N=rsqrtx(help);
  
  help=(help<1) ? sqrtx(1-help) : 0; //a1
  sincos(in[4],out[0][0]+1,out[0][0]);
  out[0][0][0]*=help;
  out[0][0][1]*=help;
  
  out[1][0][0]=in[6]; //b1
  out[1][0][1]=in[7];
  
  help=((help=1-out[0][0][0]*out[0][0][0]-out[0][0][1]*out[0][0][1]-out[1][0][0]*out[1][0][0]-out[1][0][1]*out[1][0][1])>0) ? sqrtx(help) : 0;
  sincos(in[5],out[2][0]+1,out[2][0]);
  out[2][0][0]*=help; //c1
  out[2][0][1]*=help;
  
  complexx p1={rep_N*out[2][0][0],-rep_N*out[2][0][1]}; //p1=conj(c1)/N
  complexx p2={rep_N*out[1][0][0],+rep_N*out[1][0][1]}; //p2=b1/N
  
  unsafe_complef_conj2_prod(out[1][1],out[0][1],out[0][0]); //b2
  safe_complef_prod(out[1][1],p2,out[1][1],tmp);
  complef_summ_the_conj2_prod(out[1][1],p1,out[0][2]);
  complef_prod_real(out[1][1],out[1][1],-rep_N);

  unsafe_complef_conj2_prod(out[1][2],out[0][2],out[0][0]); //b3
  safe_complef_prod(out[1][2],p2,out[1][2],tmp);
  complef_subt_the_conj2_prod(out[1][2],p1,out[0][1]);
  complef_prod_real(out[1][2],out[1][2],-rep_N);

  unsafe_complef_conj2_prod(out[2][1],out[0][1],out[0][0]); //c2
  safe_complef_conj1_prod(out[2][1],p1,out[2][1],tmp);
  complef_subt_the_conj_conj_prod(out[2][1],p2,out[0][2]);
  complef_prod_real(out[2][1],out[2][1],-rep_N);

  unsafe_complef_conj2_prod(out[2][2],out[0][2],out[0][0]); //c3
  safe_complef_conj1_prod(out[2][2],p1,out[2][2],tmp);
  complef_summ_the_conj_conj_prod(out[2][2],p2,out[0][1]);
  complef_prod_real(out[2][2],out[2][2],-rep_N);
}

