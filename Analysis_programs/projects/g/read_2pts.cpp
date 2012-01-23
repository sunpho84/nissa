#include <include.h>
#include <iostream>

using namespace std;

int nmass=2;
int T=48;
int njack=20;

int prop_combo(int r1,int im1,int r2,int im2,int ith2)
{return r1+2*(im1+nmass*(r2+2*(im2+nmass*ith2)));}

int D_combo(int r1,int r2,int th)
{return prop_combo(r1,0,r2,1,th);}

int D_combo_rev(int r1,int r2,int th)
{return prop_combo(r1,1,r2,0,th);}

jvec load_D(const char *corr,int th)
{
  int ic1=D_combo(0,0,th);
  int ic2=D_combo(1,1,th);
  int ic3=D_combo_rev(0,0,th);
  int ic4=D_combo_rev(1,1,th);
  cout<<ic1<<" "<<ic2<<endl;
  return 0.25*(jvec_load(corr,T,njack,ic1)+jvec_load(corr,T,njack,ic2)+jvec_load(corr,T,njack,ic3)+jvec_load(corr,T,njack,ic4));
}

int main()
{
  jvec eff_P5P5[2];
  for(int ith=0;ith<2;ith++)
    {
      jvec P5P5=load_D("P5P5",ith);
      eff_P5P5[ith]=effective_mass(P5P5.simmetrized(1));
      ofstream out_P5P5(combine("eff_mass_P5P5_th%d.xmg",ith).c_str());
      out_P5P5<<"@type xydy"<<endl;
      out_P5P5<<eff_P5P5[ith]<<endl;
      out_P5P5.close();
    }
  
  jvec VKVK=(load_D("V1V1",0)+load_D("V2V2",0)+load_D("V3V3",0))/3;
  jvec eff_VKVK=effective_mass(VKVK.simmetrized(1));
  ofstream out_VKVK("eff_mass_VKVK_th0.xmg");
  out_VKVK<<"@type xydy"<<endl;
  out_VKVK<<eff_VKVK<<endl;
  out_VKVK.close();
  
  
  ofstream out_diff("diff");
  out_diff<<"@type xydy"<<endl;
  out_diff<<eff_P5P5[1]-eff_VKVK<<endl;
  out_diff.close();
  
  return 0;
}
