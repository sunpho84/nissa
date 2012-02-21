#pragma once

//Refers to the doc: "doc/eo_inverter.lyx" for explenations

//apply even-odd or odd-even part of tmD, multiplied by -2
//configuration borders are assumed to have been already communicated
void tmn2Deo_or_tmn2Doe_eos(spincolor *out,spincolor *in,quad_su3 **conf,int eooe)
{
  int Xup,Xdw;
  
  if(eooe==0) communicate_od_spincolor_borders(in);
  else        communicate_ev_spincolor_borders(in);
  
  for(int X=0;X<loc_volh;X++)
    {
      color temp_c0,temp_c1,temp_c2,temp_c3;
      
      //Forward 0
      Xup=loceo_neighup[eooe][X][0];
      color_summ(temp_c0,in[Xup][0],in[Xup][2]);
      color_summ(temp_c1,in[Xup][1],in[Xup][3]);
      unsafe_su3_prod_color(out[X][0],conf[eooe][X][0],temp_c0);
      unsafe_su3_prod_color(out[X][1],conf[eooe][X][0],temp_c1);
      color_copy(out[X][2],out[X][0]);
      color_copy(out[X][3],out[X][1]);
      
      //Backward 0
      Xdw=loceo_neighdw[eooe][X][0];
      color_subt(temp_c0,in[Xdw][0],in[Xdw][2]);
      color_subt(temp_c1,in[Xdw][1],in[Xdw][3]);
      unsafe_su3_dag_prod_color(temp_c2,conf[!eooe][Xdw][0],temp_c0);
      unsafe_su3_dag_prod_color(temp_c3,conf[!eooe][Xdw][0],temp_c1);
      color_summassign(out[X][0],temp_c2);
      color_summassign(out[X][1],temp_c3);
      color_subtassign(out[X][2],temp_c2);
      color_subtassign(out[X][3],temp_c3);
      
      //Forward 1
      Xup=loceo_neighup[eooe][X][1];
      color_isumm(temp_c0,in[Xup][0],in[Xup][3]);
      color_isumm(temp_c1,in[Xup][1],in[Xup][2]);
      unsafe_su3_prod_color(temp_c2,conf[eooe][X][1],temp_c0);
      unsafe_su3_prod_color(temp_c3,conf[eooe][X][1],temp_c1);
      color_summassign(out[X][0],temp_c2);
      color_summassign(out[X][1],temp_c3);
      color_isubtassign(out[X][2],temp_c3);
      color_isubtassign(out[X][3],temp_c2);
      
      //Backward 1
      Xdw=loceo_neighdw[eooe][X][1];
      color_isubt(temp_c0,in[Xdw][0],in[Xdw][3]);
      color_isubt(temp_c1,in[Xdw][1],in[Xdw][2]);
      unsafe_su3_dag_prod_color(temp_c2,conf[!eooe][Xdw][1],temp_c0);
      unsafe_su3_dag_prod_color(temp_c3,conf[!eooe][Xdw][1],temp_c1);
      color_summassign(out[X][0],temp_c2);
      color_summassign(out[X][1],temp_c3);
      color_isummassign(out[X][2],temp_c3);
      color_isummassign(out[X][3],temp_c2);
      
      //Forward 2
      Xup=loceo_neighup[eooe][X][2];
      color_summ(temp_c0,in[Xup][0],in[Xup][3]);
      color_subt(temp_c1,in[Xup][1],in[Xup][2]);
      unsafe_su3_prod_color(temp_c2,conf[eooe][X][2],temp_c0);
      unsafe_su3_prod_color(temp_c3,conf[eooe][X][2],temp_c1);
      color_summassign(out[X][0],temp_c2);
      color_summassign(out[X][1],temp_c3);
      color_subtassign(out[X][2],temp_c3);
      color_summassign(out[X][3],temp_c2);
      
      //Backward 2
      Xdw=loceo_neighdw[eooe][X][2];
      color_subt(temp_c0,in[Xdw][0],in[Xdw][3]);
      color_summ(temp_c1,in[Xdw][1],in[Xdw][2]);
      unsafe_su3_dag_prod_color(temp_c2,conf[!eooe][Xdw][2],temp_c0);
      unsafe_su3_dag_prod_color(temp_c3,conf[!eooe][Xdw][2],temp_c1);
      color_summassign(out[X][0],temp_c2);
      color_summassign(out[X][1],temp_c3);
      color_summassign(out[X][2],temp_c3);
      color_subtassign(out[X][3],temp_c2);
      
      //Forward 3
      Xup=loceo_neighup[eooe][X][3];
      color_isumm(temp_c0,in[Xup][0],in[Xup][2]);
      color_isubt(temp_c1,in[Xup][1],in[Xup][3]);
      unsafe_su3_prod_color(temp_c2,conf[eooe][X][3],temp_c0);
      unsafe_su3_prod_color(temp_c3,conf[eooe][X][3],temp_c1);
      color_summassign(out[X][0],temp_c2);
      color_summassign(out[X][1],temp_c3);
      color_isubtassign(out[X][2],temp_c2);
      color_isummassign(out[X][3],temp_c3);
      
      //Backward 3
      Xdw=loceo_neighdw[eooe][X][3];
      color_isubt(temp_c0,in[Xdw][0],in[Xdw][2]);
      color_isumm(temp_c1,in[Xdw][1],in[Xdw][3]);
      unsafe_su3_dag_prod_color(temp_c2,conf[!eooe][Xdw][3],temp_c0);
      unsafe_su3_dag_prod_color(temp_c3,conf[!eooe][Xdw][3],temp_c1);
      color_summassign(out[X][0],temp_c2);
      color_summassign(out[X][1],temp_c3);
      color_isummassign(out[X][2],temp_c2);
      color_isubtassign(out[X][3],temp_c3);
    }
}

//wrappers
void tmn2Doe_eos(spincolor *out,spincolor *in,quad_su3 **conf){tmn2Deo_or_tmn2Doe_eos(out,in,conf,1);}
void tmn2Deo_eos(spincolor *out,spincolor *in,quad_su3 **conf){tmn2Deo_or_tmn2Doe_eos(out,in,conf,0);}

//implement ee or oo part of Dirac operator, equation(3)
void tmDee_or_oo_eos(spincolor *out,spincolor *in,double kappa,double mu)
{
  if(in==out) crash("in==out!");
  complex z={1/(2*kappa),mu};
  
  for(int X=0;X<loc_volh;X++)
    for(int ic=0;ic<3;ic++)
      {
	for(int id=0;id<2;id++) unsafe_complex_prod(out[X][id][ic],in[X][id][ic],z);
	for(int id=2;id<4;id++) unsafe_complex_conj2_prod(out[X][id][ic],in[X][id][ic],z);
      }
}

//inverse
void inv_tmDee_or_oo_eos(spincolor *out,spincolor *in,double kappa,double mu)
{
  if(in==out) crash("in==out!");
  double a=1/(2*kappa),b=mu,nrm=a*a+b*b;
  complex z={+a/nrm,-b/nrm};
  
  for(int X=0;X<loc_volh;X++)
    for(int ic=0;ic<3;ic++)
      {
	for(int id=0;id<2;id++) unsafe_complex_prod(out[X][id][ic],in[X][id][ic],z);
	for(int id=2;id<4;id++) unsafe_complex_conj2_prod(out[X][id][ic],in[X][id][ic],z);
      }
}

//implement Koo defined in equation (7) 
void tmDkern_eoprec_eos(spincolor *out,spincolor *temp,spincolor *in,quad_su3** conf,double kappa,double mu)
{
  tmn2Deo_eos(out,in,conf);
  inv_tmDee_or_oo_eos(temp,out,kappa,mu);
  tmn2Doe_eos(out,temp,conf);  
  inv_tmDee_or_oo_eos(temp,out,kappa,mu);
  tmDee_or_oo_eos(temp,in,kappa,mu);
  
  for(int ivol=0;ivol<loc_volh;ivol++)
    for(int id=0;id<2;id++)
      for(int ic=0;ic<3;ic++)
	for(int ri=0;ri<2;ri++)
	  { //gamma5 is explicitely implemented
	    out[ivol][id  ][ic][ri]=+temp[ivol][id  ][ic][ri]-out[ivol][id  ][ic][ri]*0.25;
	    out[ivol][id+2][ic][ri]=-temp[ivol][id+2][ic][ri]+out[ivol][id+2][ic][ri]*0.25;
	  }
}

//square of Koo
void tmDkern_eoprec_square_eos(spincolor *out,spincolor *temp1,spincolor *temp2,spincolor *in,quad_su3 **conf,double kappa,double mu)
{
  tmDkern_eoprec_eos(temp1,temp2, in,   conf,kappa,-mu);
  tmDkern_eoprec_eos(out,  temp2, temp1,conf,kappa,+mu);
}
