#pragma once

//Apply the dirac operator to a spincolor
//it is assumed that boundary condition has been already adjusted outside
void apply_dirac_operator(spincolor *out,spincolor *in,quad_su3 *conf,double kappac,double m)
{
  for(int ivol=0;ivol<loc_vol;ivol++)
    {
    }
}

//Apply the D+D- operator to a spincolor
//if no temporary spinor is provided, it will be allocated inside
void apply_DDdag_operator(spincolor *out,spincolor *in,spincolor *temp,quad_su3 *conf,double kappac,double m)
{
  int all=0;
  if(temp==NULL)
    {
      all=1;
      temp=(spincolor*)malloc(sizeof(spincolor)*(loc_vol+loc_bord));
    }
  
  apply_dirac_operator(temp,in,conf,kappac,+m);
  //to be added: comunication of the borders
  apply_dirac_operator(out,temp,conf,kappac,-m);
  
  if(all!=0) free(temp);
}
