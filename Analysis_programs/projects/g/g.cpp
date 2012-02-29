#include <include.h>

#define REAL 0
#define IMAG 1

int T,L,tsep;
double th;
int njack;
const double Zv[3]={0.5816,0.6103,0.6451};

int ibeta;

int tmin_V,tmax_V;
int tminL_P,tmaxL_P;
int tminS_P,tmaxS_P;
int tmin_g,tmax_g;

jvec load_3pts_charm_spec(const char *name,int reim)
{return jvec_load(combine("3pts_charm_spec_%s",name).c_str(),T,njack,reim);}

jvec load_3pts_light_spec(const char *name,int reim)
{return jvec_load(combine("3pts_light_spec_%s",name).c_str(),T,njack,reim);}

jvec load_2pts_lc(const char *name,int reim,const char *sl)
{return jvec_load(combine("2pts_lc_%s_%s",name,sl).c_str(),T,njack,reim);}

jvec load_2pts_lc_sl(const char *name,int reim)
{return load_2pts_lc(name,reim,"sl");}

jvec load_2pts_lc_ss(const char *name,int reim)
{return load_2pts_lc(name,reim,"ss");}

jack jack_average(jack &a,jack &b)
{
  double ea=a.err();
  double eb=b.err();
  
  double wa=1/(ea*ea);
  double wb=1/(eb*eb);
  
  return (a*wa+b*wb)/(wa+wb);
}

void read_input()
{
  FILE *input_file=open_file("analysis_pars","r");
  read_formatted_from_file_expecting((char*)(&L),input_file,"%d","L");
  read_formatted_from_file_expecting((char*)(&njack),input_file,"%d","njack");
  T=2*L;
  read_formatted_from_file_expecting((char*)(&ibeta),input_file,"%d","ibeta");
  read_formatted_from_file_expecting((char*)(&th),input_file,"%lg","th");
  read_formatted_from_file_expecting((char*)(&tsep),input_file,"%d","tsep");
  
  read_formatted_from_file_expecting((char*)(&tmin_V),input_file,"%d","tmin_V");
  read_formatted_from_file_expecting((char*)(&tmax_V),input_file,"%d","tmax_V");
  
  read_formatted_from_file_expecting((char*)(&tminL_P),input_file,"%d","tminL_P");
  read_formatted_from_file_expecting((char*)(&tmaxL_P),input_file,"%d","tmaxL_P");
  
  read_formatted_from_file_expecting((char*)(&tminS_P),input_file,"%d","tminS_P");
  read_formatted_from_file_expecting((char*)(&tmaxS_P),input_file,"%d","tmaxS_P");
  
  read_formatted_from_file_expecting((char*)(&tmin_g),input_file,"%d","tmin_g");
  read_formatted_from_file_expecting((char*)(&tmax_g),input_file,"%d","tmax_g");
  
  fclose(input_file);
}

int main()
{
  read_input();
  
  ///////////////////////////// Load two points for standing D and D* //////////////////////////
  
  //load ss P5P5 for D
  jvec P5P5_ss=load_2pts_lc_ss("P5P5",REAL);
  //load sl P5P5 for D
  jvec P5P5_sl=load_2pts_lc_sl("P5P5",REAL);
  
  //load ss VKVK for D
  jvec V1V1_ss=load_2pts_lc_ss("V1V1",REAL);
  jvec V2V2_ss=load_2pts_lc_ss("V2V2",REAL);
  jvec V3V3_ss=load_2pts_lc_ss("V3V3",REAL);
  jvec VKVK_ss=(V1V1_ss+V2V2_ss+V3V3_ss)/3;
  //load sl VKVK for D
  jvec V1V1_sl=load_2pts_lc_sl("V1V1",REAL);
  jvec V2V2_sl=load_2pts_lc_sl("V2V2",REAL);
  jvec V3V3_sl=load_2pts_lc_sl("V3V3",REAL);
  jvec VKVK_sl=(V1V1_sl+V2V2_sl+V3V3_sl)/3;
  
  //////////////////////////////////// Fit masses and Z for standing D and D* ////////////////////////////////////////
  
  //compute D mass and Z
  jack M_P5,ZL_P5,ZS_P5;
  two_pts_SL_fit(M_P5,ZL_P5,ZS_P5,P5P5_sl.simmetrized(1),P5P5_ss.simmetrized(1),tminL_P,tmaxL_P,tminS_P,tmaxS_P,"MSL_P5.xmg","MSS_P5.xmg");
  cout<<"D mass: "<<M_P5<<", Z: "<<ZL_P5<<endl;
  
  //compute D* mass and Z
  jack MSL_VK,ZSL_VK;
  two_pts_fit(MSL_VK,ZSL_VK,VKVK_sl.simmetrized(1),tmin_V,tmax_V,"MSL_VK.xmg","ZSL_VK.xmg");
  jack MSS_VK,ZSS_VK;
  two_pts_fit(MSS_VK,ZSS_VK,VKVK_ss.simmetrized(1),tmin_V,tmax_V,"MSS_VK.xmg","ZSS_VK.xmg");
  jack M_VK=jack_weighted_average(MSS_VK,MSL_VK),ZS_VK=sqrt(ZSS_VK),ZL_VK=ZSL_VK/ZS_VK;
  
  //reconstuct moving D mass
  double qi=M_PI*th/L;
  double q2=3*sqr(qi);
  jack Eth_P5=sqrt(sqr(M_P5)+q2);
  
  //reconstruct numerically and semi-analitically the time dependance of three points
  jvec Dth_DV_td_sa(T,njack),Dth_DV_td_nu(T,njack);
  for(int t=0;t<=tsep;t++)
    {
      Dth_DV_td_nu[t]=P5P5_sl[(tsep-t)%T]*VKVK_sl[t]/(ZL_P5*ZL_VK);
      Dth_DV_td_sa[t]=(ZS_P5*ZS_VK)*exp(-(tsep-t)*Eth_P5)*exp(-(t)*M_VK)/(2*Eth_P5*2*M_VK);
    }
  for(int t=tsep+1;t<T;t++)
    {
      Dth_DV_td_sa[t]=P5P5_sl[(t-tsep)%T]*VKVK_sl[T-t]/(ZL_P5*ZL_VK);
      Dth_DV_td_nu[t]=(ZS_P5*ZS_VK)*exp(-(t-tsep)*Eth_P5)*exp(-(T-t)*M_VK)/(2*Eth_P5*2*M_VK);
    }
  
  //compare time dependance and its simmetric
  if(tsep==T/2)
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
  
  ///////////////////////////// Determine the matrix element of AK between D(th=+-) and D* ////////////////////////////
  
  jvec Dth_AK_DV_sa=P5thAKVK/Dth_DV_td_sa;
  jvec Dth_AK_DV_nu=P5thAKVK/Dth_DV_td_nu;
  
  //compare matrix element and its simmetric
  if(tsep==T/2)
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
  jack R1_sa=constant_fit( (tsep==T/2) ? Dth_AK_DV_sa.simmetrized(1) : Dth_AK_DV_sa,tmin_g,tmax_g,"R1_sa.xmg");
  jack R1_nu=constant_fit( (tsep==T/2) ? Dth_AK_DV_nu.simmetrized(1) : Dth_AK_DV_nu,tmin_g,tmax_g,"R1_nu.xmg");
  
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
  if(tsep==T/2)
    {
      ofstream out("R2_pt1_simm.xmg");
      out<<"@type xydy"<<endl;
      out<<R2_pt1_corr<<endl;
      
      out<<"@type xydy"<<endl;
      out<<R2_pt1_corr.simmetrized(1)<<endl;
    }
  
  //write part 2 of R2 simmetrized and not
  if(tsep==T/2)
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
  jack gc_pt1=fpi_gDvDPi_pt1*Zv[ibeta]/(2*sqrt(M_P5*M_VK));
  jack gc_pt2=gc_pt1*fpi_gDvDPi_corr_rel;
  
  //print gc
  cout<<"gc: \n pt1: "<<gc_pt1<<"\n pt2: "<<gc_pt2<<endl;
  
  
  
  return 0;
}
