#include "nissa.hpp"

using namespace nissa;

int main()
{
  double d;
  while(fread(&d,sizeof(double),1,stdin)==1)
    {
      change_endianness(d);
      if(fwrite(&d,sizeof(double),1,stdout)!=1) CRASH("writing to stdout");
    }
  
  return 0;
}
