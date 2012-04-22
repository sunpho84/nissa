#include "nissa.h"

int main()
{
  check_endianess();
  printf("Endianess of this machine: ");
  if(little_endian==1) printf("little\n");
  else printf("big\n");
  
  return 0;
}
