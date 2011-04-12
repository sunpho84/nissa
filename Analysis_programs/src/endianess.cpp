#pragma once

//check the endianess of the machine
int check_endianess()
{
  int big_endian=1;
  big_endian=(int)(*(char*)(&big_endian));
  
  return big_endian;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////

//fucking tool to revert the endianess of doubles
void doubles_to_doubles_changing_endianess(double *dest,double *sour,int ndoubles)
{
  char *cdest,*csour;
  char temp;
  
  if(dest==sour)
    for(int idouble=0;idouble<ndoubles;idouble++)
      {
        cdest=(char*)(dest+idouble);
        csour=(char*)(sour+idouble);
        
        temp=csour[7];
        csour[7]=cdest[0];
        cdest[0]=csour[7];

        temp=csour[6];
        csour[6]=cdest[1];
        cdest[1]=csour[6];

        temp=csour[5];
        csour[5]=cdest[2];
        cdest[2]=csour[5];

        temp=csour[4];
        csour[4]=cdest[3];
        cdest[3]=csour[4];
      }
  else
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
