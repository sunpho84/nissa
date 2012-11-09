#ifndef _DIRAC_OPERATOR_WSTAT_H
#define _DIRAC_OPERATOR_WSTAT_H
void apply_Wstat(spincolor *out,quad_su3 *conf,spincolor *in,int,int);
void apply_Wstat2(spincolor *out,quad_su3 *conf,spincolor *temp,spincolor *in);
void reconstruct_Wstat(spincolor *outminus,spincolor *outplus,quad_su3 *conf,spincolor *in);
#endif
