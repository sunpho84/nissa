#include <common.cpp>

int main()
{
  define_combo_info();
  
  jvec V0P5=log(-load_corr_3pts("V0P5",STAND_DSOUR_POSTH_DSEQ_CSPEC,REAL).simmetrized(-1));
  jvec test(T/2,njack);
  for(int t=0;t<T/2;t++)
    test[t]=V0P5[t]-V0P5[t+1];

  ofstream out("Dth_minus_D.xmg");
  out<<"@type xydy"<<endl;
  out<<test<<endl;
  return 0;
}
