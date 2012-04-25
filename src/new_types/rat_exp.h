#ifndef _RAT_EXP_H
#define _RAT_EXP_H

void master_printf_rat_approx(rat_approx *appr);
void rat_approx_create(rat_approx *appr,int nterms,const char *name);
void rat_approx_destroy(rat_approx *appr);

#endif
