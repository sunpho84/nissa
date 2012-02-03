#include <common.cpp>

void compute_ZV(info_3pts &three,info_2pts &two)
{
  jvec V0P5=load_corr_3pts("V0P5",three,REAL).simmetrized(-1);
  jvec P5P5=load_corr_2pts("P5P5",two,REAL);
  
  jack cenP5P5=P5P5[L]/2;
  
  {
    ofstream out("V0P5.xmg");
    out<<"@type xydy"<<endl;
    out<<V0P5<<endl;
  }
  
  {
    ofstream out("P5P5.xmg");
    out<<"@type xydy"<<endl;
    out<<effective_mass(P5P5.simmetrized(1))<<endl;
  }
  
  jack ZV=-cenP5P5/constant_fit(V0P5,8,16);
  cout<<ZV<<endl;
}

int main()
{
  define_combo_info();
  
  compute_ZV(POSTH_D_POSTH_D_CSPEC,POSTH_D_CSPEC_00);
  compute_ZV(STAND_D_STAND_D_CSPEC,STAND_D_CSPEC_00);
  compute_ZV(STAND_ETA_STAND_ETA,STAND_ETA_00);
  compute_ZV(STAND_D_STAND_D_LSPEC,STAND_D_LSPEC_00);
  
  return 0;
}
