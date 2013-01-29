#include <string.h>
#include <stdlib.h>

#include "new_types_definitions.h"
#include "../base/routines.h"

//allocate a new rat approx
void rat_approx_create(rat_approx_type *appr,int nterms,const char *name)
{
  memcpy(appr->name,name,20);
  appr->minimum=appr->maximum=0;
  appr->nterms=nterms;
  appr->poles=(double*)malloc(sizeof(double)*2*nterms);
  appr->weights=appr->poles+nterms;
}

//free a rational approx
void rat_approx_destroy(rat_approx_type *appr)
{
  free(appr->poles);
}

//print a rational approximation
void master_printf_rat_approx(rat_approx_type *appr)
{
  master_printf("Rational approximation %s of x^%lg:\n",appr->name,appr->exp_power);
  master_printf("  valid in the interval: %lg %lg\n",appr->minimum,appr->maximum);
  master_printf("  const: %lg\n",appr->cons);
  master_printf("  nterms: %d\n",appr->nterms);
  for(int iterm=0;iterm<appr->nterms;iterm++)
    master_printf("   %d) pole: %lg, weight: %lg\n",iterm,appr->poles[iterm],appr->weights[iterm]);
}
