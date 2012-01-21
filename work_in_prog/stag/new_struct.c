#pragma once

typedef struct
{
  int deg;
  double mass;
  double re_pot;
  double im_pot;
  double charge;
} quark_content;

//allocate a new rat approx
void rat_approx_create(rat_approx *appr,int nterms,char *name)
{
  memcpy(appr->name,name,20);
  appr->nterms=nterms;
  appr->poles=malloc(sizeof(double)*2*nterms);
  appr->weights=appr->poles+nterms;
}

//free a rational approx
void rat_approx_destroy(rat_approx *appr)
{
  free(appr->poles);
}

//print a rational approximation
void master_printf_rat_approx(rat_approx *appr)
{
  master_printf("Rational approximation %s of x^%lg:\n",appr->name,appr->exp_power);
  master_printf("  valid in the interval: %lg %lg\n",appr->minimum,appr->maximum);
  master_printf("  const: %lg\n",appr->cons);
  master_printf("  nterms: %d\n",appr->nterms);
  for(int iterm=0;iterm<appr->nterms;iterm++)
    master_printf("   %d) pole: %lg, weight: %lg\n",iterm,appr->poles[iterm],appr->weights[iterm]);
}
