#include <include.h>
#include <iostream>

using namespace std;

int main()
{
  jtvec c(jvec_load("/home/francesco/QCD/LAVORI/RHO/sim2/000/data",16,10,8*2+2));
  jtvec s(c.simmetrized(1));
  jtvec f(effective_mass(s));
  
  cout<<c<<endl<<s<<endl<<f;
  
  jack E=constant_fit(f,3,6),Z;
  
  cout<<E<<endl;
  
  linear_fit(f,E,Z,3,6);
  cout<<E<<" "<<Z<<endl;
  
  return 0;
}
