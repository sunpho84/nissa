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
double color_norm2_diff(color *a,color *b)
{
  double loc_norm2=0;
  double norm2;
  
  nissa_loc_volh_loop(ivol)
    for(int ic=0;ic<3;ic++)
      for(int ri=0;ri<2;ri++)
	{	
	  double d=a[ivol][ic][ri]-b[ivol][ic][ri];
	  loc_norm2+=d*d;
	}
  MPI_Allreduce(&loc_norm2,&norm2,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  
  return norm2;
}
