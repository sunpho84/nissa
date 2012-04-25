#ifndef _DIRAC_OPERATOR_TMQ2_H
#define _DIRAC_OPERATOR_TMQ2_H

void apply_tmQ2(spincolor *out,quad_su3 *conf,double kappa,spincolor *temp,double mu,spincolor *in);
void apply_tmQ2_RL(spincolor *out,quad_su3 *conf,double kappa,spincolor *temp,int RL,double mu,spincolor *in);
void apply_tmQ2_left(spincolor *out,quad_su3 *conf,double kappa,spincolor *temp,double mu,spincolor *in);
void apply_tmQ2_m2_RL(spincolor *out,quad_su3 *conf,double kappa,spincolor *temp,int RL,double m2,spincolor *in);

#endif
