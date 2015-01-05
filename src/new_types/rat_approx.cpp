#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <string.h>
#include <stdlib.h>

#include "new_types_definitions.hpp"

#include "routines/ios.hpp"

namespace nissa
{
  //allocate a new rat approx
  void rat_approx_create(rat_approx_t *appr,int degree,const char *name)
  {
    memcpy(appr->name,name,20);
    appr->minimum=appr->maximum=0;
    appr->degree=degree;
    appr->poles=(double*)malloc(sizeof(double)*2*degree);
    appr->weights=appr->poles+degree;
  }
  
  //free a rational approx
  void rat_approx_destroy(rat_approx_t *appr)
  {
    free(appr->poles);
  }
  
  //print a rational approximation
  void master_printf_rat_approx(rat_approx_t *appr)
  {
    master_printf("Rational approximation %s of x^(%d/%d):\n",appr->name,appr->num,appr->den);
    master_printf("  valid in the interval: %lg %lg\n",appr->minimum,appr->maximum);
    master_printf("  const: %lg\n",appr->cons);
    master_printf("  degree: %d\n",appr->degree);
    for(int iterm=0;iterm<appr->degree;iterm++)
      master_printf("   %d) pole: %lg, weight: %lg\n",iterm,appr->poles[iterm],appr->weights[iterm]);
  }
}
