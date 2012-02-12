#include <include.h>

#define REAL 0
#define IMAG 1

const int T=48,L=T/2;
const int njack=16;

int tmin_V=10,tmax_V=15;
int tmin_P=7,tmax_P=23;
int tmin_g=tmin_P,tmax_g=min(T/2-tmin_P,tmax_V);

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
  jvec P5thAKVJ=-(P5thA1V2+P5thA1V3+P5thA2V1+P5thA2V3+P5thA3V1+P5thA3V2)/6;
  
  //put together the equation (16)
  jvec P5thAKVJK=P5thAKVK+P5thAKVJ;

  //compare the different parts
  {
    ofstream out("P5thAKVJK_parts.xmg");
    out<<"@type xydy"<<endl;
    out<<P5thAKVK<<endl;
    out<<"&/n@type xydy"<<endl;
    out<<P5thAKVJ<<endl;
    out<<"&/n@type xydy"<<endl;
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
  jack Mth_P5=sqrt(sqr(M_P5)+q2);
  
  //reconstruct semi-analitically the time dependance of three points
  jvec Dth_DV_td(T,njack);
  for(int t=0;t<T/2;t++)
    {
      Dth_DV_td[t    ]=sqrt(Z2_P5*Z2_VK)*exp(-t*M_VK)*exp(-(T/2-t)*Mth_P5)/(2*Mth_P5*2*M_VK);
      Dth_DV_td[t+T/2]=sqrt(Z2_P5*Z2_VK)*exp(-t*Mth_P5)*exp(-(T/2-t)*M_VK)/(2*Mth_P5*2*M_VK);
    }
  
  ///////////////////////////// Determine the matrix element of AK between D(th=+-) and D* ////////////////////////////
  
  jvec Dth_AK_DV=P5thAKVK/Dth_DV_td;
  
  //compare matrix element and its simmetric
  {
    ofstream out("Dth_AK_DV_simm.xmg");
    out<<"@type xydy"<<endl;
    out<<Dth_AK_DV<<endl;
    out<<"&/n@type xydy"<<endl;
    out<<Dth_AK_DV.simmetric()<<endl;
  }
  
  
  /////////////////////////////// Compute g ///////////////////////////
  
  //determine g correlation function
  jvec g_corr=(Dth_AK_DV*0.746/(M_P5+M_VK)).simmetrized(1);
  
  //fit g
  jack g=constant_fit(g_corr,tmin_g,tmax_g);
  
  //write matrix element
  {
    ofstream out("g_corr_numerical.xmg");
    out<<"@type xydy"<<endl;
    out<<g_corr<<endl;
    out<<"&"<<endl<<"@type xy"<<endl;
    out<<tmin_g<<" "<<g.med()+g.err()<<endl;
    out<<tmax_g<<" "<<g.med()+g.err()<<endl;
    out<<tmax_g<<" "<<g.med()-g.err()<<endl;
    out<<tmin_g<<" "<<g.med()-g.err()<<endl;
    out<<tmin_g<<" "<<g.med()+g.err()<<endl;
  }
  
  //print g
  cout<<"g: "<<g<<endl;
  
  
  /////////////////////////////// Load the three points for the corrections /////////////////////////////////////

  int TEST=IMAG;
  jvec P5thA0V1_00=load_3pts_charm_spec("A0V1",TEST);
  jvec P5thA0V2_00=load_3pts_charm_spec("A0V2",TEST);
  jvec P5thA0V3_00=load_3pts_charm_spec("A0V3",TEST);
  jvec P5thA0VK_00=(P5thA0V1_00+P5thA0V2_00+P5thA0V3_00)/3;
  
  //write
  {
    ofstream out("P5thA0VK.xmg");
    out<<"@type xydy"<<endl;
    out<<P5thA0VK_00/P5thAKVJK<<endl;
    
    out<<"@type xydy"<<endl;
    out<<(P5thA0VK_00/P5thAKVJK).simmetrized(1)<<endl;
  }  
  
  cout<<(M_VK*M_VK-M_P5*M_P5)/(2*qi*M_VK)<<endl;
  cout<<M_VK<<endl<<M_P5<<endl;
  
  return 0;
}
