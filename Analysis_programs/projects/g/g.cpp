#include <include.h>

#define REAL 0
#define IMAG 1

const int T=48,L=T/2;
const int njack=16;
const double ZV=0.746;

int tmin_V=7,tmax_V=15;
int tmin_P=7,tmax_P=23;
int tmin_g=12,tmax_g=15;

jvec load_3pts_charm_spec(const char *name,int reim)
{return jvec_load(combine("3pts_charm_spec_%s",name).c_str(),T,njack,reim);}

jvec load_3pts_light_spec(const char *name,int reim)
{return jvec_load(combine("3pts_light_spec_%s",name).c_str(),T,njack,reim);}

jvec load_2pts_lc(const char *name,int reim)
{return jvec_load(combine("2pts_lc_%s_ss",name).c_str(),T,njack,reim);}


int main()
{
  //////////////////////////////// Calculate C2 defined in eq. (16) ////////////////////////////////

  //load C2 for D(th=+-) and D*, with spectator c
  //this is done in two parts, both for positive and negative th
  
  //load 1st part
  jvec P5thA1V1=load_3pts_charm_spec("A1V1",REAL);
  jvec P5thA2V2=load_3pts_charm_spec("A2V2",REAL);
  jvec P5thA3V3=load_3pts_charm_spec("A3V3",REAL);
  jvec P5thAKVK=(P5thA1V1+P5thA2V2+P5thA3V3)/3;

  //load 2nd part
  jvec P5thA1V2=load_3pts_charm_spec("A1V2",REAL);
  jvec P5thA1V3=load_3pts_charm_spec("A1V3",REAL);
  jvec P5thA2V1=load_3pts_charm_spec("A2V1",REAL);
  jvec P5thA2V3=load_3pts_charm_spec("A2V3",REAL);
  jvec P5thA3V1=load_3pts_charm_spec("A3V1",REAL);
  jvec P5thA3V2=load_3pts_charm_spec("A3V2",REAL);
  jvec P5thAKVJ=(P5thA1V2+P5thA1V3+P5thA2V1+P5thA2V3+P5thA3V1+P5thA3V2)/6;
  
  //put together the equation (16)
  jvec P5thAKVJK=P5thAKVK-P5thAKVJ;

  //compare the different parts
  {
    ofstream out("P5thAKVJK_parts.xmg");
    out<<"@type xydy"<<endl;
    out<<P5thAKVK<<endl;
    out<<"&\n@type xydy"<<endl;
    out<<-P5thAKVJ<<endl;
    out<<"&\n@type xydy"<<endl;
    out<<P5thAKVJK<<endl;
  }
  
  ///////////////////////////// Load two points for standing D and D* //////////////////////////
  
  //load P5P5 for D
  jvec P5P5=load_2pts_lc("P5P5",REAL);
  
  //load VKVK for D
  jvec V1V1=load_2pts_lc("V1V1",REAL);
  jvec V2V2=load_2pts_lc("V2V2",REAL);
  jvec V3V3=load_2pts_lc("V3V3",REAL);
  jvec VKVK=(V1V1+V2V2+V3V3)/3;
  
  //////////////////////////////////// Fit masses and Z for standing D and D* ////////////////////////////////////////
  
  //compute D mass and Z
  jack M_P5,Z2_P5;
  P5P5_fit(M_P5,Z2_P5,P5P5.simmetrized(1),tmin_P,tmax_P,"M_P5.xmg","Z2_P5.xmg");
  
  //compute D* mass and Z
  jack M_VK,Z2_VK;
  P5P5_fit(M_VK,Z2_VK,VKVK.simmetrized(1),tmin_V,tmax_V,"M_VK.xmg","Z2_VK.xmg");
  
  //reconstuct moving D mass
  double th=0.41;
  double qi=M_PI*th/L;
  double q2=3*sqr(qi);
  jack Eth_P5=sqrt(sqr(M_P5)+q2);
  
  //reconstruct numerically and semi-analitically the time dependance of three points
  jvec Dth_DV_td_sa(T,njack),Dth_DV_td_nu(T,njack);
  for(int t=0;t<=T/2;t++)
    {
      Dth_DV_td_nu[t    ]=P5P5[(T/2+t)%T]*VKVK[t]/sqrt(Z2_P5*Z2_VK);
      Dth_DV_td_sa[t    ]=sqrt(Z2_P5*Z2_VK)*exp(-t*M_VK)*exp(-(T/2-t)*Eth_P5)/(2*Eth_P5*2*M_VK);
    }
  for(int t=1;t<T/2;t++)
    {
      Dth_DV_td_sa[T/2+t]=Dth_DV_td_sa[T/2-t];
      Dth_DV_td_nu[T/2+t]=Dth_DV_td_nu[T/2-t];
    }
  
  //compare time dependance and its simmetric
  {
    ofstream out("Dth_DV_td_sa_simm.xmg");
    out<<"@type xydy"<<endl;
    out<<Dth_DV_td_sa<<endl;
    out<<"&\n@type xydy"<<endl;
    out<<Dth_DV_td_sa.simmetric()<<endl;
  }
  
  
  //compare time dependance semi-analytical and numeric
  {
    ofstream out("Dth_DV_td_sa_nu.xmg");
    out<<"@type xydy"<<endl;
    out<<Dth_DV_td_sa<<endl;
    out<<"&\n@type xydy"<<endl;
    out<<Dth_DV_td_nu<<endl;
  }
  
  
  cout<<Dth_DV_td_sa[0][njack]<<" "<<Dth_DV_td_sa[0]<<endl;
  cout<<Dth_DV_td_sa[47][njack]<<endl;
  ///////////////////////////// Determine the matrix element of AK between D(th=+-) and D* ////////////////////////////
  
  jvec Dth_AK_DV_sa=P5thAKVK/Dth_DV_td_sa;
  jvec Dth_AK_DV_nu=P5thAKVK/Dth_DV_td_nu;
  
  //compare matrix element and its simmetric
  {
    ofstream out("Dth_AK_DV_sa_simm.xmg");
    out<<"@type xydy"<<endl;
    out<<Dth_AK_DV_sa<<endl;
    out<<"&\n@type xydy"<<endl;
    out<<Dth_AK_DV_sa.simmetric()<<endl;
    out<<"&\n@type xydy"<<endl;
    out<<Dth_AK_DV_sa.simmetrized(1)<<endl;
  }
  
  
  //compare matrix element semi-analytical and numeric
  {
    ofstream out("Dth_AK_DV_sa_nu.xmg");
    out<<"@type xydy"<<endl;
    out<<Dth_AK_DV_sa<<endl;
    out<<"&\n@type xydy"<<endl;
    out<<Dth_AK_DV_nu<<endl;
  }
  
  //fit matrix element
  jack R1_sa=constant_fit(Dth_AK_DV_sa.simmetrized(1),tmin_g,tmax_g,"R1_sa.xmg");
  jack R1_nu=constant_fit(Dth_AK_DV_nu.simmetrized(1),tmin_g,tmax_g,"R1_nu.xmg");
  
  //determine the form factor
  jack A1=R1_sa/(M_VK+M_P5);
  
  /////////////////////////////// Load the three points for the corrections /////////////////////////////////////

  int TEST=REAL;
  jvec P5thA0V1=load_3pts_charm_spec("A0V1",TEST);
  jvec P5thA0V2=load_3pts_charm_spec("A0V2",TEST);
  jvec P5thA0V3=load_3pts_charm_spec("A0V3",TEST);
  jvec P5thA0VK=(P5thA0V1+P5thA0V2+P5thA0V3)/3;
  
  
  //build the ratio
  jvec R2_pt1_corr=P5thA0VK/P5thAKVJK;
  jvec R2_pt2_corr=P5thAKVJ/P5thAKVJK;

  //write part 1 of R2 simmetrized and not
  {
    ofstream out("R2_pt1_simm.xmg");
    out<<"@type xydy"<<endl;
    out<<R2_pt1_corr<<endl;
    
    out<<"@type xydy"<<endl;
    out<<R2_pt1_corr.simmetrized(1)<<endl;
  }
  
  //write part 2 of R2 simmetrized and not
  {
    ofstream out("R2_pt2_simm.xmg");
    out<<"@type xydy"<<endl;
    out<<R2_pt2_corr<<endl;
    
    out<<"@type xydy"<<endl;
    out<<R2_pt2_corr.simmetrized(1)<<endl;
  }
  
  //determine the correction
  jack R2_pt1=constant_fit(R2_pt1_corr.simmetrized(1),tmin_g,tmax_g,"R2_pt1.xmg");
  jack R2_pt2=constant_fit(R2_pt2_corr.simmetrized(1),tmin_g,tmax_g,"R2_pt2.xmg");
  jack R2_pt2_coef=(Eth_P5-M_VK)/qi;
  cout<<"R2: \n pt1: "<<R2_pt1<<"\n pt2: "<<R2_pt2<<"\n coef2: "<<R2_pt2_coef<<endl;
  jack R2=-(R2_pt1+R2_pt2*R2_pt2_coef);
  
  jack A2frA1=R2*sqr(M_VK+M_P5)/(2*qi*M_VK);
  cout<<"A2/A1: "<<A2frA1<<endl;
  
  //////////////// Compute the full value of fpi*gDvDPi ////////////////
  
  jack fpi_gDvDPi_pt1=(M_VK+M_P5)*A1;
  jack fpi_gDvDPi_corr_rel=(M_VK-M_P5)/(M_VK+M_P5)*A2frA1;

  
  /////////////////////////////// Compute gc ///////////////////////////
  
  //determine gc
  jack gc_pt1=fpi_gDvDPi_pt1*ZV/(2*sqrt(M_P5*M_VK));
  jack gc_pt2=gc_pt1*fpi_gDvDPi_corr_rel;
  
  //print gc
  cout<<"gc: \n pt1: "<<gc_pt1<<"\n pt2: "<<gc_pt2<<endl;
  
  
  
  return 0;
}
