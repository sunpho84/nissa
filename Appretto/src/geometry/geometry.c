#pragma once

#include "geometry_lx.c"
#include "geometry_eo.c"

//swap the data from lexical order to even odd
void swap_lx_to_eo_or_eo_to_lx(char *vect_e,char *vect_o,char *vect_lx,int nbytes_per_site,int bord,int eotolx_lxtoeo)
{
  int tot_site=loc_vol;
  char *vect_eo[2]={vect_e,vect_o};

  if(bord) tot_site+=loc_bord;

  for(int loclx=0;loclx<tot_site;loclx++)
    {
      int loceo=loceo_of_loclx[loclx];
      int par=loclx_parity[loclx];
      
      if(eotolx_lxtoeo) memcpy(vect_lx+nbytes_per_site*loclx,vect_eo[par]+nbytes_per_site*loceo,nbytes_per_site);
      else              memcpy(vect_eo[par]+nbytes_per_site*loceo,vect_lx+nbytes_per_site*loclx,nbytes_per_site);
    }
}

//swap the data from even odd to lexical order
void swap_lx_to_eo(char *out_e,char *out_o,char *in_lx,int nbytes_per_site,int bord)
{
  swap_lx_to_eo_or_eo_to_lx(out_e,out_o,in_lx,nbytes_per_site,bord,0);
}
void swap_eo_to_lx(char *out_lx,char *in_e,char *in_o,int nbytes_per_site,int bord)
{
  swap_lx_to_eo_or_eo_to_lx(in_e,in_o,out_lx,nbytes_per_site,bord,1);
}

//wrappers
void swap_spincolor_lx_to_eo(spincolor *out_e,spincolor *out_o,spincolor *in_lx,int bord)
{
  swap_lx_to_eo((char*)out_e,(char*)out_o,(char*)in_lx,sizeof(spincolor),bord);
}
void swap_spincolor_eo_to_lx(spincolor *out_lx,spincolor *in_e,spincolor *in_o,int bord)
{
  swap_eo_to_lx((char*)out_lx,(char*)in_e,(char*)in_o,sizeof(spincolor),bord);
}
