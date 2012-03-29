#include <stdio.h>
#include <math.h>
#include "quad.c"

int main()
{
  double du=1,ds=1;
  for(int i=1;i<73;i++)
    if(i%2==0) ds+=du/(2*i+1);
    else ds-=du/(2*i+1);

  printf("%.18lg\n",4*ds);
  
  __float128 qu=sqrt(2),qs=sqrt(2);
    
  printf("%.18lg\n",(double)(qu*qu-2));
  
  return 0;
}
