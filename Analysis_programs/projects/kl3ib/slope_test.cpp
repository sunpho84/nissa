#include "include.h"
#include "kl3_common.cpp"

jvec ratio_P5P5_fit(jvec ratio,jvec corr)
{
  int njack=corr.njack;
  int nel=corr.nel,TH=nel-1;
  
  //fit the mass
  jack M,Z;
  P5P5_fit(M,Z,corr,12,23,"M_P5P5","Z");
  cout<<M<<" "<<Z<<endl;
  
  //fit the slope
  jvec diff(nel-1,njack);
  for(int t=0;t<nel-1;t++)
    for(int ij=0;ij<njack+1;ij++)
      {
	diff[t].data[ij]=ratio[t+1].data[ij]-ratio[t].data[ij];
	diff[t].data[ij]/=(t+1-TH)*tanh(M[ij]*(t+1-TH))-(t-TH)*tanh(M[ij]*(t-TH));
      }
  jack SL=constant_fit(diff,12+1,23,"slope");
  cout<<SL<<endl;
  
  //fit the constant
  diff=jvec(nel,njack);
  for(int t=0;t<nel;t++)
    for(int ij=0;ij<njack+1;ij++)
      diff[t].data[ij]=ratio[t].data[ij]-SL[ij]*(t-TH)*tanh(M[ij]*(t-TH));
  jack A=constant_fit(diff,12,23,"Constant");
  cout<<A<<endl;
  
  return diff;
}

jvec load_chaveraged_2pts(const char *nome_file,int im1,int im2,int ik1,int ik2,int ri)
{
  char path[1024];
  sprintf(path,"%s/%s",base_path,nome_file);  
  
  //read the 2 two points with opposite r
  return (read_two_points(path,im1,im2,ik1,ik2,0,0,ri)+
	  read_two_points(path,im1,im2,ik1,ik2,1,1,ri))/2;
}

int main()
{
  njack=10;
  read_input();
  
  int im1=0,im2=1; //quello raddoppiato e' il primo!
  int ik1=0,ik2=0;
  int ri=0;
  
  //load the data
  jvec Pi_P5P5 =load_chaveraged_2pts("oPPo-ss_conf.1.dat",im1,im1,ik1,ik2,ri).simmetrized(1);
  jvec K_A0P5 =load_chaveraged_2pts("oA0Po-ss_conf.1.dat",im1,im2,ik1,ik2,ri).simmetrized(-1);
  jvec Ks_A0P5=load_chaveraged_2pts("oA0Po-sd_conf.1.dat",im1,im2,ik1,ik2,ri).simmetrized(-1);
  jvec K_P5P5 =load_chaveraged_2pts("oPPo-ss_conf.1.dat",im1,im2,ik1,ik2,ri).simmetrized(1);
  jvec Ks_P5P5=load_chaveraged_2pts("oPPo-sd_conf.1.dat",im1,im2,ik1,ik2,ri).simmetrized(1);
  
  //define the ratios
  jvec ratio_P5P5=Ks_P5P5/K_P5P5;
  jvec ratio_A0P5=Ks_A0P5/K_A0P5;
 
  ratio_P5P5_fit(ratio_P5P5,K_P5P5);
  
  return 0;
}
