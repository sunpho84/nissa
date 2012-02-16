#include <common.cpp>

int tmin_V=10,tmax_V=15;
int tmin_P=7,tmax_P=23;
int tmin_g=tmin_P,tmax_g=min(T/2-tmin_P,tmax_V);

int main()
{
  define_combo_info();

  //////////////////////////////// Calculate C2 defined in eq. (16) ////////////////////////////////

  //load C2 for D(th=+-) and D*, with spectator c
  //this is done in two parts, both for positive and negative th
  
  //load 1st part for positive th
  jvec P5pA1V1=load_corr_3pts("A1V1",STAND_DSOUR_POSTH_DSEQ_CSPEC,REAL);
  jvec P5pA2V2=load_corr_3pts("A2V2",STAND_DSOUR_POSTH_DSEQ_CSPEC,REAL);
  jvec P5pA3V3=load_corr_3pts("A3V3",STAND_DSOUR_POSTH_DSEQ_CSPEC,REAL);
  jvec P5pAKVK=(P5pA1V1+P5pA2V2+P5pA3V3)/3;

  //load 1st part for negative th
  jvec P5mA1V1=load_corr_3pts("A1V1",STAND_DSOUR_POSTH_DSEQ_CSPEC,REAL);
  jvec P5mA2V2=load_corr_3pts("A2V2",STAND_DSOUR_NEGTH_DSEQ_CSPEC,REAL);
  jvec P5mA3V3=load_corr_3pts("A3V3",STAND_DSOUR_NEGTH_DSEQ_CSPEC,REAL);
  jvec P5mAKVK=(P5mA1V1+P5mA2V2+P5mA3V3)/3;
  
  //load 2nd part for positive th
  jvec P5pA1V2=load_corr_3pts("A1V2",STAND_DSOUR_POSTH_DSEQ_CSPEC,REAL);
  jvec P5pA1V3=load_corr_3pts("A1V3",STAND_DSOUR_POSTH_DSEQ_CSPEC,REAL);
  jvec P5pA2V1=load_corr_3pts("A2V1",STAND_DSOUR_POSTH_DSEQ_CSPEC,REAL);
  jvec P5pA2V3=load_corr_3pts("A2V3",STAND_DSOUR_POSTH_DSEQ_CSPEC,REAL);
  jvec P5pA3V1=load_corr_3pts("A3V1",STAND_DSOUR_POSTH_DSEQ_CSPEC,REAL);
  jvec P5pA3V2=load_corr_3pts("A3V2",STAND_DSOUR_POSTH_DSEQ_CSPEC,REAL);
  jvec P5pAKVJ=-(P5pA1V2+P5pA1V3+P5pA2V1+P5pA2V3+P5pA3V1+P5pA3V2)/6;
  
  //load 2nd part for positive th
  jvec P5mA1V2=load_corr_3pts("A1V2",STAND_DSOUR_NEGTH_DSEQ_CSPEC,REAL);
  jvec P5mA1V3=load_corr_3pts("A1V3",STAND_DSOUR_NEGTH_DSEQ_CSPEC,REAL);
  jvec P5mA2V1=load_corr_3pts("A2V1",STAND_DSOUR_NEGTH_DSEQ_CSPEC,REAL);
  jvec P5mA2V3=load_corr_3pts("A2V3",STAND_DSOUR_NEGTH_DSEQ_CSPEC,REAL);
  jvec P5mA3V1=load_corr_3pts("A3V1",STAND_DSOUR_NEGTH_DSEQ_CSPEC,REAL);
  jvec P5mA3V2=load_corr_3pts("A3V2",STAND_DSOUR_NEGTH_DSEQ_CSPEC,REAL);
  jvec P5mAKVJ=-(P5mA1V2+P5mA1V3+P5mA2V1+P5mA2V3+P5mA3V1+P5mA3V2)/6;
  
  //put together the equation (16) for positive and negative th
  jvec P5pAKVJK=P5pAKVK+P5pAKVJ;
  jvec P5mAKVJK=P5mAKVK+P5mAKVJ;

  //compare the different parts for positive th
  {
    ofstream out("P5pAKVJK_parts.xmg");
    out<<"@type xydy"<<endl;
    out<<P5pAKVK<<endl;
    out<<"&/n@type xydy"<<endl;
    out<<P5pAKVJ<<endl;
    out<<"&/n@type xydy"<<endl;
    out<<P5pAKVJK<<endl;
  }
  
  //compare the different parts for negative th
  {
    ofstream out("P5mAV_parts.xmg");
    out<<"@type xydy"<<endl;
    out<<P5mAKVK<<endl;
    out<<"&/n@type xydy"<<endl;
    out<<P5mAKVJ<<endl;
    out<<"&/n@type xydy"<<endl;
    out<<P5mAKVJK<<endl;
  }
  
  ///////////////////////////// Load two points for standing D and D*, and moving D //////////////////////////
  
  //load P5P5 for D with spectator c
  jvec P5P5_00=load_corr_2pts("P5P5",STAND_D_CSPEC_00,REAL);
  jvec P5P5_11=load_corr_2pts("P5P5",STAND_D_CSPEC_11,REAL);
  jvec P5P5=(P5P5_00+P5P5_11)/2;
  
  //load VKVK for D with spectator c
  jvec V1V1_00=load_corr_2pts("V1V1",STAND_D_CSPEC_00,REAL);
  jvec V2V2_00=load_corr_2pts("V2V2",STAND_D_CSPEC_00,REAL);
  jvec V3V3_00=load_corr_2pts("V3V3",STAND_D_CSPEC_00,REAL);
  jvec VKVK_00=(V1V1_00+V2V2_00+V3V3_00)/3;
  jvec V1V1_11=load_corr_2pts("V1V1",STAND_D_CSPEC_11,REAL);
  jvec V2V2_11=load_corr_2pts("V2V2",STAND_D_CSPEC_11,REAL);
  jvec V3V3_11=load_corr_2pts("V3V3",STAND_D_CSPEC_11,REAL);
  jvec VKVK_11=(V1V1_11+V2V2_11+V3V3_11)/3;
  jvec VKVK=(VKVK_00+VKVK_11)/2;
  
  //load moving D
  jvec P5P5p_00=load_corr_2pts("P5P5",POSTH_D_CSPEC_00,REAL);
  jvec P5P5m_00=load_corr_2pts("P5P5",NEGTH_D_CSPEC_00,REAL);
  jvec P5P5pm_00=(P5P5p_00+P5P5m_00)/2;
  
  
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
  jack Mpm_P5=sqrt(sqr(M_P5)+q2);
  
  //compute moving D mass and Z
  jack Mpm_P5_num,Z2pm_P5;
  P5P5_fit(Mpm_P5_num,Z2pm_P5,P5P5pm_00.simmetrized(1),tmin_P,tmax_P,"Mpm_P5_num.xmg","Z2pm_P5.xmg");  

  //reconstruct semi-analitically the time dependance of three points
  jvec D_DV_td(T,njack);
  for(int t=0;t<T/2;t++)
    {
      D_DV_td[t    ]=sqrt(Z2_P5*Z2_VK)*exp(-t*M_VK)*exp(-(T/2-t)*Mpm_P5)/(2*Mpm_P5*2*M_VK);
      D_DV_td[t+T/2]=sqrt(Z2_P5*Z2_VK)*exp(-t*Mpm_P5)*exp(-(T/2-t)*M_VK)/(2*Mpm_P5*2*M_VK);
    }
  //reconstruct numerically the time dependance of three points
  jvec D_DV_td_num(T,njack);
  for(int t=0;t<T/2;t++)
    {
      D_DV_td_num[t    ]=P5P5pm_00[T/2-t]*VKVK[t    ]/sqrt(Z2pm_P5*Z2_VK);
      D_DV_td_num[t+T/2]=P5P5pm_00[t    ]*VKVK[T/2-t]/sqrt(Z2pm_P5*Z2_VK);
    }
  
  //write the t dependance both numerical and semi-analitical
  {
    ofstream out("D_DV_td.xmg");
    out<<"@type xydy"<<endl;
    out<<D_DV_td<<endl;
    out<<"&/n@type xydy"<<endl;
    out<<D_DV_td_num<<endl;
  }

  
  ///////////////////////////// Determine the matrix element of AK between D(th=+-) and D* ////////////////////////////
  
  jvec Dp_AK_DV=P5pAKVK/D_DV_td;
  jvec Dm_AK_DV=P5mAKVK/D_DV_td;
  
  //average positive and negative th
  jvec Dpm_AK_DV=(Dp_AK_DV+Dm_AK_DV)/2;
  
  //compare positive and negative th matrix element with their average
  {
    ofstream out("D_AK_DV.xmg");
    out<<"@type xydy"<<endl;
    out<<Dp_AK_DV<<endl;
    out<<"&/n@type xydy"<<endl;
    out<<Dm_AK_DV<<endl;
    out<<"&/n@type xydy"<<endl;
    out<<Dpm_AK_DV<<endl;
  }
  
  //compare matrix element and its simmetric
  {
    ofstream out("D_AK_DV_simm.xmg");
    out<<"@type xydy"<<endl;
    out<<Dpm_AK_DV<<endl;
    out<<"&/n@type xydy"<<endl;
    out<<Dpm_AK_DV.simmetric()<<endl;
  }
  
  
  //////////////////////////////////// Compute g ///////////////////////////
  
  //load VKVK for D with spectator c
  //determine g correlation function
  jvec g_corr=(Dpm_AK_DV*0.746/(M_P5+M_VK)).simmetrized(1);
  
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
  jvec P5pA0V1_00=load_corr_3pts("A0V1",STAND_DSOUR_POSTH_DSEQ_CSPEC,TEST);
  jvec P5pA0V2_00=load_corr_3pts("A0V2",STAND_DSOUR_POSTH_DSEQ_CSPEC,TEST);
  jvec P5pA0V3_00=load_corr_3pts("A0V3",STAND_DSOUR_POSTH_DSEQ_CSPEC,TEST);
  jvec P5pA0VK_00=(P5pA0V1_00+P5pA0V2_00+P5pA0V3_00)/3;
  jvec P5mA0V1_00=load_corr_3pts("A0V1",STAND_DSOUR_NEGTH_DSEQ_CSPEC,TEST);
  jvec P5mA0V2_00=load_corr_3pts("A0V2",STAND_DSOUR_NEGTH_DSEQ_CSPEC,TEST);
  jvec P5mA0V3_00=load_corr_3pts("A0V3",STAND_DSOUR_NEGTH_DSEQ_CSPEC,TEST);
  jvec P5mA0VK_00=(P5mA0V1_00+P5mA0V2_00+P5mA0V3_00)/3;
  
  //write
  {
    ofstream out("P5pA0VK.xmg");
    out<<"@type xydy"<<endl;
    out<<P5pA0VK_00/P5pAKVJK<<endl;
  }  
  //write
  {
    ofstream out("P5mA0VK.xmg");
    out<<"@type xydy"<<endl;
    out<<P5mA0VK_00/P5mAKVJK<<endl;
  }
  //write
  {
    ofstream out("P5A0VK.xmg");
    out<<"@type xydy"<<endl;
    
    jvec ciao=(((P5pA0VK_00/P5pAKVJK)+(P5mA0VK_00/P5mAKVJK))/2).simmetrized(1);
    cout<<constant_fit(ciao,tmin_g,tmax_g)<<endl;
    out<<ciao<<endl;
  }  

  cout<<(M_VK*M_VK-M_P5*M_P5)/(2*qi*M_VK)<<endl;
  cout<<M_VK<<endl<<M_P5<<endl;
  return 0;
}
