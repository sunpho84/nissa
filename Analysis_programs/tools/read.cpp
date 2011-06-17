#include <include.h>
#include <iostream>

using namespace std;

int main()
{
  int T=48;
  int njack=16;
  
  jvec a=jvec_load("/home/francesco/QCD/LAVORI/NF2/DATA_LOW/3.90/24/0.0064/VKVK",T,njack,17).simmetrized(1);
  jvec b=jvec_load("/home/francesco/QCD/LAVORI/NF2/DATA/3.90/24/0.0064/VKVK",T,njack,28).simmetrized(1);
  jvec c=(a+b)/2;
  jvec ea=effective_mass(a);
  jvec eb=effective_mass(b);
  jvec ec=effective_mass(c);
  jvec ed=(ea+eb)/2;
  
  ofstream out("/tmp/out");
  out<<"@type xydy"<<endl;
  out<<ea<<endl;
  out<<"&"<<endl;
  out<<eb<<endl;
  out<<"&"<<endl;
  out<<ec<<endl;

  jack ma=constant_fit(ea,11,23);
  jack mb=constant_fit(eb,11,23);
  jack mc=constant_fit(ec,11,23);
  jack md=constant_fit(ed,11,23);
  jack me=(ma+mb)/2;
  
  cout<<ma<<endl;
  cout<<mb<<endl;
  cout<<mc<<endl;
  cout<<md<<endl;
  cout<<me<<endl;
  
  return 0;
}
