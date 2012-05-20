#include <string.h>

#include "../base/global_variables.h"
#include "../base/vectors.h"

//swap the data from lexical order to even odd
void swap_lx_to_eo_or_eo_to_lx(char *vect_e,char *vect_o,char *vect_lx,int nbytes_per_site,int bord,int eotolx_lxtoeo)
{
  int tot_site=loc_vol;
  char *vect_eo[2]={vect_e,vect_o};

  if(bord) tot_site+=bord_vol;

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
{swap_lx_to_eo_or_eo_to_lx(out_e,out_o,in_lx,nbytes_per_site,bord,0);}
void swap_eo_to_lx(char *out_lx,char *in_e,char *in_o,int nbytes_per_site,int bord)
{swap_lx_to_eo_or_eo_to_lx(in_e,in_o,out_lx,nbytes_per_site,bord,1);}

//wrappers
void swap_spincolor_lx_to_eo(spincolor *out_e,spincolor *out_o,spincolor *in_lx,int bord)
{swap_lx_to_eo((char*)out_e,(char*)out_o,(char*)in_lx,sizeof(spincolor),bord);}
void swap_spincolor_eo_to_lx(spincolor *out_lx,spincolor *in_e,spincolor *in_o,int bord)
{swap_eo_to_lx((char*)out_lx,(char*)in_e,(char*)in_o,sizeof(spincolor),bord);}

//separate the even or odd part of a vector
void take_e_or_o_part_of_lx_vector(char *out_e_or_o,char *in_lx,int bps,int par)
{
  //extract
  nissa_loc_vol_loop(loclx)
    if(par==loclx_parity[loclx])
      memcpy(out_e_or_o+bps*loceo_of_loclx[loclx],in_lx+bps*loclx,bps);
  
  set_borders_invalid(out_e_or_o);
}
//wrappers
void take_e_part_of_lx_vector(char *out_e,char *in_lx,int bps)
{take_e_or_o_part_of_lx_vector(out_e,in_lx,bps,EVN);}
void take_o_part_of_lx_vector(char *out_o,char *in_lx,int bps)
{take_e_or_o_part_of_lx_vector(out_o,in_lx,bps,ODD);}
void take_e_part_of_lx_color(color *out_e,color *in_lx)
{take_e_part_of_lx_vector((char*)out_e,(char*)in_lx,sizeof(color));}
void take_o_part_of_lx_color(color *out_o,color *in_lx)
{take_o_part_of_lx_vector((char*)out_o,(char*)in_lx,sizeof(color));}

//separate the even and odd part of a vector
void split_lx_vector_into_eo_parts(char **out_eo,char *in_lx,int bps)
{
  //split
  nissa_loc_vol_loop(loclx)
    memcpy(out_eo[loclx_parity[loclx]]+bps*loceo_of_loclx[loclx],in_lx+bps*loclx,bps);
  
  set_borders_invalid(out_eo[0]);
  set_borders_invalid(out_eo[1]);
}

//paste the even and odd parts of a vector into a full lx vector
void paste_eo_parts_into_lx_vector(char *out_lx,char **in_eo,int bps)
{
  //paste
  nissa_loc_vol_loop(loclx)
    {
      int eo=loceo_of_loclx[loclx];
      int par=loclx_parity[loclx];
      memcpy(out_lx+bps*loclx,in_eo[par]+bps*eo,bps);
    }

  set_borders_invalid(out_lx);
}

//wrappers
void split_lx_conf_into_eo_parts(quad_su3 **eo_out,quad_su3 *lx_in)
{split_lx_vector_into_eo_parts((char**)eo_out,(char*)lx_in,sizeof(quad_su3));}
void split_lx_color_into_eo_parts(color **eo_out,color *lx_in)
{split_lx_vector_into_eo_parts((char**)eo_out,(char*)lx_in,sizeof(color));}
void split_lx_spincolor_into_eo_parts(spincolor **eo_out,spincolor *lx_in)
{split_lx_vector_into_eo_parts((char**)eo_out,(char*)lx_in,sizeof(spincolor));}
void split_lx_spin_into_eo_parts(spin **eo_out,spin *lx_in)
{split_lx_vector_into_eo_parts((char**)eo_out,(char*)lx_in,sizeof(spin));}

void paste_eo_parts_into_lx_conf(quad_su3 *out_lx,quad_su3 **in_eo)
{paste_eo_parts_into_lx_vector((char*)out_lx,(char**)in_eo,sizeof(quad_su3));}
void paste_eo_parts_into_lx_color(color *out_lx,color **in_eo)
{paste_eo_parts_into_lx_vector((char*)out_lx,(char**)in_eo,sizeof(color));}
void paste_eo_parts_into_lx_spincolor(spincolor *out_lx,spincolor **in_eo)
{paste_eo_parts_into_lx_vector((char*)out_lx,(char**)in_eo,sizeof(spincolor));}
void paste_eo_parts_into_lx_spin(spin *out_lx,spin **in_eo)
{paste_eo_parts_into_lx_vector((char*)out_lx,(char**)in_eo,sizeof(spin));}
