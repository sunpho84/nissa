#pragma once

//fucking tool to revert the endianess of doubles
void revert_endianess_double_vector(double *dest,double *sour,int ndoubles)
{
  char *cdest,*csour;

  for(int idouble=0;idouble<ndoubles;idouble++)
    {
      cdest=(char*)(dest+idouble);
      csour=(char*)(sour+idouble);

      cdest[0]=csour[7];
      cdest[1]=csour[6];
      cdest[2]=csour[5];
      cdest[3]=csour[4];
      cdest[4]=csour[3];
      cdest[5]=csour[2];
      cdest[6]=csour[1];
      cdest[7]=csour[0];
    }
}
