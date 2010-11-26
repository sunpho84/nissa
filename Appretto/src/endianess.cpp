#pragma once

#include <iostream>

using namespace std;

//check the endianess of the machine
void check_endianess()
{
  big_endian=1;
  big_endian=int(*(char*)&big_endian);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////

//fucking tool to revert the endianess of doubles
void doubles_to_doubles_changing_endianess(double *dest,double *sour,int ndoubles)
{
  char *cdest,*csour;
  
  if(rank==0 and debug>1) cout<<"Reverting the endianess ot the data"<<endl;
  
  for(int idouble=0;idouble<ndoubles;idouble++)
    {
      cdest=(char*)(dest+idouble);
      csour=(char*)(sour+idouble);
      
      swap(cdest[0],csour[7]);
      swap(cdest[1],csour[6]);
      swap(cdest[2],csour[5]);
      swap(cdest[3],csour[4]);
    }
}

////////////////////Copy a vector of floats to doubles. Sweep is reversed to avoid overwriting////////////////

//Do not change endianess
void floats_to_doubles_same_endianess(double *dest,float *sour,int n)
{
  if(rank==0 and debug>1) cout<<"Converting "<<n<<" floats to doubles"<<endl;
  
  for(int i=n-1;i>=0;i--) dest[i]=(double)sour[i];
}

//Change endianess
void floats_to_doubles_changing_endianess(double *dest,float *sour,int n)
{
  char *c;

  if(rank==0 and debug>1) cout<<"Converting "<<n<<" floats to doubles changing endianess"<<endl;

  for(int i=n-1;i>=0;i--)
    {
      c=(char*)(sour+i);
      
      swap(c[0],c[3]);
      swap(c[1],c[2]);
      
      dest[i]=(double)sour[i];
    }
}

////////////////////Copy a vector of doubles to floats. Sweep is direct, to avoid overwriting////////////////

//Do not change the endianess
void doubles_to_floats_same_endianess(float *dest,double *sour,int n)
{
  if(rank==0 and debug>1) cout<<"Converting "<<n<<" doubles to floats"<<endl;

  for(int i=0;i<n;i++) dest[i]=(float)sour[i];
}

//Change endianess
void doubles_to_floats_changing_endianess(float *dest,double *sour,int n)
{
  char *c;

  if(rank==0 and debug>1) cout<<"Converting "<<n<<" doubles to floats changing endianess"<<endl;

  for(int i=0;i<n;i++)
    {
      dest[i]=(double)sour[i];

      c=(char*)(sour+i);

      swap(c[0],c[3]);
      swap(c[1],c[2]);
    }
}
