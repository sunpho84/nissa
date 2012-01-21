#pragma once

//debug
void read_e_color(color *e,char *path)
{
  color *lx=nissa_malloc("lx",loc_vol,color);
  
  read_color(lx,path);
  take_e_part_of_lx_color(e,lx,loc_vol);
  
  nissa_free(lx);
}
void read_o_color(color *o,char *path)
{
  color *lx=nissa_malloc("lx",loc_vol,color);
  
  read_color(lx,path);
  take_o_part_of_lx_color(o,lx,loc_vol);
  
  nissa_free(lx);
}
