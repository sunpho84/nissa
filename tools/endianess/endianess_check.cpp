#include "nissa.hpp"

using namespace nissa;

int main()
{
  check_endianness();
  printf("Endianess of this machine: ");
  if(little_endian==1) printf("little\n");
  else printf("big\n");
  
  return 0;
}
