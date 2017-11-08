#include "nissa.hpp"

using namespace nissa;

int main()
{
  double d;
  while(fread(&d,sizeof(double),1,stdin)==1)
    {
      change_endianness(d);
      fwrite(&d,sizeof(double),1,stdin);
    }
  
  return 0;
}
