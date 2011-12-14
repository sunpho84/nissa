#include "nissa.h"

int main()
{
  check_endianess();
  printf("Endianess of this machine: ");
  if(big_endian==1) printf("big\n");
  else printf("little\n");
  
  return 0;
}
