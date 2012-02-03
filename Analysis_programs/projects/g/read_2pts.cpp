#include <include.h>
#include <iostream>

using namespace std;

int nspec=2;
int ntheta=3;
int nmass=2;
int T=48,L=24;
int njack=13;

int prop_combo(int r1,int im1,int ispec,int r2,int im2,int ith2,int reim)
{return reim+2*(r1+2*(im1+nmass*(r2+2*(im2+nmass*(ith2+ntheta*ispec)))));}

int D_combo(int r1,int ispec,int r2,int th,int reim)
{return prop_combo(r1,0,ispec,r2,1,th,reim);}

int D_combo_rev(int r1,int ispec,int r2,int th,int reim)
{return prop_combo(r1,1,ispec,r2,0,th,reim);}

jvec load_D(const char *corr,int th,int reim)
{
  int ic1=D_combo(0,0,0,th,reim);
  int ic2=D_combo(1,0,1,th,reim);
  int ic3=D_combo_rev(0,0,0,th,reim);
  int ic4=D_combo_rev(1,0,1,th,reim);
  cout<<ic1<<" "<<ic2<<endl;
  return 0.25*(jvec_load(corr,T,njack,ic1)+jvec_load(corr,T,njack,ic2)+jvec_load(corr,T,njack,ic3)+jvec_load(corr,T,njack,ic4));
}

int main()
{
  jvec eff_P5P5[3];
  for(int ith=0;ith<3;ith++)
    {
      jvec P5P5=load_D("2pts_P5P5_30_30",ith,0);
      eff_P5P5[ith]=effective_mass(P5P5.simmetrized(1));
      ofstream out_P5P5(combine("eff_mass_P5P5_th%d.xmg",ith).c_str());
      out_P5P5<<"@type xydy"<<endl;
      out_P5P5<<eff_P5P5[ith]<<endl;
      out_P5P5.close();
    }
  
  jvec VKVK=(load_D("2pts_V1V1_30_30",0,0)+load_D("2pts_V2V2_30_30",0,0)+load_D("2pts_V3V3_30_30",0,0))/3;
  jvec eff_VKVK=effective_mass(VKVK.simmetrized(1));
  ofstream out_VKVK("eff_mass_VKVK_th0.xmg");
  out_VKVK<<"@type xydy"<<endl;
  out_VKVK<<eff_VKVK<<endl;
  out_VKVK.close();
  
  double q=M_PI*0.5/L*sqrt(3);
  ofstream out_diff("diff");
  out_diff<<"@type xydy"<<endl;
  out_diff<<sqrt(sqr(eff_VKVK)-sqr(eff_P5P5[1]+eff_P5P5[2])/4)-q<<endl;
  out_diff.close();
  
  ofstream out_test("test");
  out_test<<"@type xydy"<<endl;
  out_test<<(eff_P5P5[1]+eff_P5P5[2])/2-eff_P5P5[0]<<endl;
  out_test<<"&\n@type xydy"<<endl;
  out_test<<sqr(eff_VKVK-eff_P5P5[0])/(2*eff_VKVK)<<endl;
  out_test.close();
  
  return 0;
}
