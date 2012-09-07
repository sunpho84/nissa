#ifndef _DIRAC_OPERATOR_STD_H
#define _DIRAC_OPERATOR_STD_H

void apply_st2Doe(color *out,quad_su3 **conf,color *in);
void apply_stD2ee(color *out,quad_su3 **conf,color *temp,double mass,color *in);
void apply_stDeo_half(color *out,quad_su3 **conf,color *in);
void apply_stDoe(color *out,quad_su3 **conf,color *in);
void apply_stD2ee_m2(color *out,quad_su3 **conf,color *temp,double m2,color *in);
void apply_stD2ee_zero_mass(color *out,quad_su3 **conf,color *temp,color *in);
void evn_apply_stD(color *out,quad_su3 **conf,double m,color **in);

#endif
