#pragma once

#include <mpi.h>
#include <stdio.h>

#include "geometry_lx.c"

//swap the data from lexical order to even odd
void unsafe_swap_lx_to_eo(char *out,char *in,int nbytes_per_site)
{
  for(int i=0;i<loc_vol;i++)
    {
      memcpy(out+(i/2+(i%2)*loc_volr)*nbytes_per_site,in,nbytes_per_site);
      in+=nbytes_per_site;
    }
}

//swap the data from even odd to lexical order
void unsafe_swap_eo_to_lx(char *out,char *in,int nbytes_per_site)
{
  for(int i=0;i<loc_vol;i++)
    {
      memcpy(out+(i/(loc_volr)+(i%loc_volr)*2)*nbytes_per_site,in,nbytes_per_site);
      in+=nbytes_per_site;
    }
}

//wrappers
void unsafe_swap_gauge_lx_to_eo(quad_su3 *out,quad_su3 *in)
{
  unsafe_swap_lx_to_eo((char*)out,(char*)in,sizeof(quad_su3));
}

void unsafe_swap_gauge_eo_to_lx(quad_su3 *out,quad_su3 *in)
{
  unsafe_swap_eo_to_lx((char*)out,(char*)in,sizeof(quad_su3));
}

