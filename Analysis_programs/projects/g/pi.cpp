#include <include.h>

#define REAL 0
#define IMAG 1

int T,L,tsep;
double th;
const int njack=16;

int ibeta;

int tmin_P,tmax_P;
int tmin_V,tmax_V;

jvec load_2pts_ll(const char *name,int reim,const char *sl)
{return jvec_load(combine("2pts_ll_%s_%s",name,sl).c_str(),T,njack,reim);}

jvec load_2pts_ll_sl(const char *name,int reim)
{return load_2pts_ll(name,reim,"sl");}

jvec load_2pts_ll_ss(const char *name,int reim)
{return load_2pts_ll(name,reim,"ss");}

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
  T=2*L;
  read_formatted_from_file_expecting((char*)(&ibeta),input_file,"%d","ibeta");
  read_formatted_from_file_expecting((char*)(&th),input_file,"%lg","th");
  read_formatted_from_file_expecting((char*)(&tsep),input_file,"%d","tsep");
  
  read_formatted_from_file_expecting((char*)(&tmin_V),input_file,"%d","tmin_V");
  read_formatted_from_file_expecting((char*)(&tmax_V),input_file,"%d","tmax_V");
  
  read_formatted_from_file_expecting((char*)(&tmin_P),input_file,"%d","tmin_P");
  read_formatted_from_file_expecting((char*)(&tmax_P),input_file,"%d","tmax_P");
  
  fclose(input_file);
}

int main()
{
  read_input();
  
  ///////////////////////////// Load two points for standing Pion //////////////////////////
  
  //load ss P5P5 for Pi
  jvec P5P5_ss=load_2pts_ll_ss("P5P5",REAL);
  //load sl P5P5 for Pi
  jvec P5P5_sl=load_2pts_ll_sl("P5P5",REAL);
  
  //load ss VKVK for rho
  jvec V1V1_ss=load_2pts_ll_ss("V1V1",REAL);
  jvec V2V2_ss=load_2pts_ll_ss("V2V2",REAL);
  jvec V3V3_ss=load_2pts_ll_ss("V3V3",REAL);
  jvec VKVK_ss=(V1V1_ss+V2V2_ss+V3V3_ss)/3;
  //load sl VKVK for rho
  jvec V1V1_sl=load_2pts_ll_sl("V1V1",REAL);
  jvec V2V2_sl=load_2pts_ll_sl("V2V2",REAL);
  jvec V3V3_sl=load_2pts_ll_sl("V3V3",REAL);
  jvec VKVK_sl=(V1V1_sl+V2V2_sl+V3V3_sl)/3;
  
  //////////////////////////////////// Fit masses and Z for standing pion and rho ////////////////////////////////////////
  
  //compute Pi mass and Z
  jack MSS_P5,ZSS_P5;
  P5P5_fit(MSS_P5,ZSS_P5,P5P5_ss.simmetrized(1),tmin_P,tmax_P,"MSS_P5.xmg","ZSS_P5.xmg");
  jack MSL_P5,ZSL_P5;
  P5P5_fit(MSL_P5,ZSL_P5,P5P5_sl.simmetrized(1),tmin_P,tmax_P,"MSL_P5.xmg","ZSL_P5.xmg");
  jack M_P5=jack_average(MSS_P5,MSL_P5),ZS_P5=sqrt(ZSS_P5),ZL_P5=ZSL_P5/ZS_P5;
  
  //compute Pi mass and Z
  jack MSS_VK,ZSS_VK;
  P5P5_fit(MSS_VK,ZSS_VK,VKVK_ss.simmetrized(1),tmin_V,tmax_V,"MSS_VK.xmg","ZSS_VK.xmg");
  jack MSL_VK,ZSL_VK;
  P5P5_fit(MSL_VK,ZSL_VK,VKVK_sl.simmetrized(1),tmin_V,tmax_V,"MSL_VK.xmg","ZSL_VK.xmg");
  jack M_VK=jack_average(MSS_VK,MSL_VK),ZS_VK=sqrt(ZSS_VK),ZL_VK=ZSL_VK/ZS_VK;
  
  return 0;
}
