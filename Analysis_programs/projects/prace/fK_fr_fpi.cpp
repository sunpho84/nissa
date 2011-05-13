#include <include.h>
#include <iostream>

#include "prace_common.cpp"

using namespace std;

int tmin[3]={10,13,15};
int q[10]={0,1,1,1,2,2,2,2,2,2};

jack calculate_f(int iml,int imq,const char *tag_l)
{
  int ru=0;
  
  //load
  jvec cDql=(read_ch_thimpr_P5_P5(ith_spec,iml,ru,imq)+read_ch_thimpr_P5_P5(ith_spec,iml,!ru,imq))/2;
  
  //fit E
  jvec temp=effective_mass(cDql.simmetrized(1)),temp2(T,njack);
  temp.print_to_file("corr_Dq%s_%d",tag_l,imq);
  jack E=constant_fit(temp,tmin[q[imq]],21);
  
  //fit Z
  for(int t=0;t<T;t++) temp2[t]=cDql[t]/exp(-E*TH)/cosh(E*(TH-t))*E;
  temp2=temp2.simmetrized(1);
  temp2.print_to_file("Z_Dq%s_%d",tag_l,imq);
  jack Z=constant_fit(temp2,tmin[q[imq]],21);
  
  return (mass[iml]+mass[imq])*sqrt(Z)/(E*sinh(E));
}

int main()
{
  read_data_list();
  
  int iml=0,ims=3;
  
  for(int imq=0;imq<nmass;imq++)
    {
      jack fDl=calculate_f(iml,imq,"l");
      jack fDs=calculate_f(ims,imq,"s");
      
      cout<<mass[imq]<<" "<<fDs/fDl<<endl;
    }
  
  return 0;
}
