#include <common.cpp>

int tmin=7,tmax=14;

int main()
{
  define_combo_info();

  //////////////////////////////// Calculate C2 defined in eq. (16) ////////////////////////////////

  //load P5AKVK for D(th=+-) and D*, with spectator c
  //this is done in two parts, both for positive and negative th
  
  //positive th
  int TEST=REAL;//IMAG;
  jvec P5pV1V2_cspec=load_corr_3pts("V1V2",STAND_DSOUR_POSTH_DSEQ_CSPEC,TEST);
  jvec P5pV2V3_cspec=load_corr_3pts("V2V3",STAND_DSOUR_POSTH_DSEQ_CSPEC,TEST);
  jvec P5pV3V1_cspec=load_corr_3pts("V3V1",STAND_DSOUR_POSTH_DSEQ_CSPEC,TEST);
  jvec P5pV2V1_cspec=load_corr_3pts("V2V1",STAND_DSOUR_POSTH_DSEQ_CSPEC,TEST);
  jvec P5pV3V2_cspec=load_corr_3pts("V3V2",STAND_DSOUR_POSTH_DSEQ_CSPEC,TEST);
  jvec P5pV1V3_cspec=load_corr_3pts("V1V3",STAND_DSOUR_POSTH_DSEQ_CSPEC,TEST);
  jvec Cpij_cspec=(P5pV1V2_cspec+P5pV2V3_cspec+P5pV3V1_cspec-P5pV2V1_cspec-P5pV3V2_cspec-P5pV1V3_cspec)/6;
  
  jvec P5mV1V2_cspec=load_corr_3pts("V1V2",STAND_DSOUR_NEGTH_DSEQ_CSPEC,TEST);
  jvec P5mV2V3_cspec=load_corr_3pts("V2V3",STAND_DSOUR_NEGTH_DSEQ_CSPEC,TEST);
  jvec P5mV3V1_cspec=load_corr_3pts("V3V1",STAND_DSOUR_NEGTH_DSEQ_CSPEC,TEST);
  jvec P5mV2V1_cspec=load_corr_3pts("V2V1",STAND_DSOUR_NEGTH_DSEQ_CSPEC,TEST);
  jvec P5mV3V2_cspec=load_corr_3pts("V3V2",STAND_DSOUR_NEGTH_DSEQ_CSPEC,TEST);
  jvec P5mV1V3_cspec=load_corr_3pts("V1V3",STAND_DSOUR_NEGTH_DSEQ_CSPEC,TEST);
  jvec Cmij_cspec=(P5mV1V2_cspec+P5mV2V3_cspec+P5mV3V1_cspec-P5mV2V1_cspec-P5mV3V2_cspec-P5mV1V3_cspec)/6;
  
  jvec P5pV1V2_lspec=load_corr_3pts("V1V2",STAND_DSOUR_POSTH_DSEQ_LSPEC,TEST);
  jvec P5pV2V3_lspec=load_corr_3pts("V2V3",STAND_DSOUR_POSTH_DSEQ_LSPEC,TEST);
  jvec P5pV3V1_lspec=load_corr_3pts("V3V1",STAND_DSOUR_POSTH_DSEQ_LSPEC,TEST);
  jvec P5pV2V1_lspec=load_corr_3pts("V2V1",STAND_DSOUR_POSTH_DSEQ_LSPEC,TEST);
  jvec P5pV3V2_lspec=load_corr_3pts("V3V2",STAND_DSOUR_POSTH_DSEQ_LSPEC,TEST);
  jvec P5pV1V3_lspec=load_corr_3pts("V1V3",STAND_DSOUR_POSTH_DSEQ_LSPEC,TEST);
  jvec Cpij_lspec=(P5pV1V2_lspec+P5pV2V3_lspec+P5pV3V1_lspec-P5pV2V1_lspec-P5pV3V2_lspec-P5pV1V3_lspec)/6;
  
  jvec P5mV1V2_lspec=load_corr_3pts("V1V2",STAND_DSOUR_NEGTH_DSEQ_LSPEC,TEST);
  jvec P5mV2V3_lspec=load_corr_3pts("V2V3",STAND_DSOUR_NEGTH_DSEQ_LSPEC,TEST);
  jvec P5mV3V1_lspec=load_corr_3pts("V3V1",STAND_DSOUR_NEGTH_DSEQ_LSPEC,TEST);
  jvec P5mV2V1_lspec=load_corr_3pts("V2V1",STAND_DSOUR_NEGTH_DSEQ_LSPEC,TEST);
  jvec P5mV3V2_lspec=load_corr_3pts("V3V2",STAND_DSOUR_NEGTH_DSEQ_LSPEC,TEST);
  jvec P5mV1V3_lspec=load_corr_3pts("V1V3",STAND_DSOUR_NEGTH_DSEQ_LSPEC,TEST);
  jvec Cmij_lspec=(P5mV1V2_lspec+P5mV2V3_lspec+P5mV3V1_lspec-P5mV2V1_lspec-P5mV3V2_lspec-P5mV1V3_lspec)/6;
  
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
  
  //compute D mass and Z
  jack M_P5,Z2_P5;
  P5P5_fit(M_P5,Z2_P5,P5P5.simmetrized(1),tmin,tmax,"M_P5.xmg","Z2_P5.xmg");
  
  //compute D* mass and Z
  jack M_VK,Z2_VK;
  P5P5_fit(M_VK,Z2_VK,VKVK.simmetrized(1),tmin,tmax,"M_VK.xmg","Z2_VK.xmg");
  
  //reconstuct moving D mass
  double th=0.41;
  jack Mpm_P5=sqrt(sqr(M_P5)+3*sqr(M_PI*th/L));

  //reconstruct semi-analitically the time dependance of three points
  jvec D_DV_td(T,njack);
  for(int t=0;t<T/2;t++)
    {
      D_DV_td[t    ]=sqrt(Z2_P5*Z2_VK)*exp(-t*M_VK)*exp(-(T/2-t)*Mpm_P5)/(2*Mpm_P5*2*M_VK);
      D_DV_td[t+T/2]=sqrt(Z2_P5*Z2_VK)*exp(-t*Mpm_P5)*exp(-(T/2-t)*M_VK)/(2*Mpm_P5*2*M_VK);
    }
  
  //write the t dependance both numerical and semi-analitical
  {
    ofstream out("D_DV_td.xmg");
    out<<"@type xydy"<<endl;
    out<<D_DV_td<<endl;
  }

  //compare the different parts for positive th
  {
    ofstream out("Cij.xmg");
    out<<"@type xydy"<<endl;
    out<<(Cpij_cspec+Cmij_cspec)/2/D_DV_td<<endl;
    //out<<"\n@type xydy"<<endl;
    //out<<Cpij_lspec<<endl;
    //out<<"\n@type xydy"<<endl;
    //out<<Cmij_lspec<<endl;
    //out<<"\n@type xydy"<<endl;
    //out<<(Cpij_lspec+Cmij_lspec)/2<<endl;
  }
  
  return 0;
}
