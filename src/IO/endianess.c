#pragma once

//check the endianess of the machine
void check_endianess()
{
  little_endian=1;
  little_endian=(int)(*(char*)(&little_endian));
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////

//fucking tool to revert the endianess of doubles
void doubles_to_doubles_changing_endianess(double *dest,double *sour,int ndoubles)
{
  char *cdest,*csour;
  char temp;
  
  if(debug_lvl>1) master_printf("Reverting the endianess ot the data\n");
  
  if(dest==sour)
    for(int idouble=0;idouble<ndoubles;idouble++)
      {
	cdest=(char*)(dest+idouble);
	csour=(char*)(sour+idouble);
	
	temp=csour[7];
	csour[7]=cdest[0];
	cdest[0]=temp;

	temp=csour[6];
	csour[6]=cdest[1];
	cdest[1]=temp;

	temp=csour[5];
	csour[5]=cdest[2];
	cdest[2]=temp;

	temp=csour[4];
	csour[4]=cdest[3];
	cdest[3]=temp;
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

void floats_to_floats_changing_endianess(float *dest,float *sour,int nint)
{
  char *cdest,*csour;
  char temp;
  
  if(debug_lvl>1) master_printf("Reverting the endianess ot the data\n");

  if(dest==sour)
    for(int iint=0;iint<nint;iint++)
      {
	cdest=(char*)(dest+iint);
	csour=(char*)(sour+iint);
	
	temp=csour[3];
	csour[3]=cdest[0];
	cdest[0]=temp;

	temp=csour[2];
	csour[2]=cdest[1];
	cdest[1]=temp;
    }
  else
    for(int iint=0;iint<nint;iint++)
      {
	cdest=(char*)(dest+iint);
	csour=(char*)(sour+iint);
	
	cdest[0]=csour[3];
	cdest[1]=csour[2];
	cdest[2]=csour[1];
	cdest[3]=csour[0];
    }
}

////////////////////Copy a vector of floats to doubles. Sweep is reversed to avoid overwriting////////////////

//Do not change endianess
void floats_to_doubles_same_endianess(double *dest,float *sour,int n)
{
  if(rank==0 && debug_lvl>1) printf("Converting %d floats to doubles\n",n);
  
  for(int i=n-1;i>=0;i--) dest[i]=(double)(sour[i]);
}

//Change endianess
void floats_to_doubles_changing_endianess(double *dest,float *sour,int n)
{
  char *c;
  char temp;

  if(rank==0 && debug_lvl>1) printf("Converting %d floats to doubles changing endianess\n",n);

  for(int i=n-1;i>=0;i--)
    {
      c=(char*)(sour+i);
      
      temp=c[3];
      c[3]=c[0];
      c[0]=temp;

      temp=c[2];
      c[2]=c[1];
      c[1]=temp;
      
      dest[i]=(double)(sour[i]);
    }
}

////////////////////Copy a vector of doubles to floats. Sweep is direct, to avoid overwriting////////////////

//Do not change the endianess
void doubles_to_floats_same_endianess(float *dest,double *sour,int n)
{
  if(rank==0 && debug_lvl>1) printf("Converting %d doubles to floats\n",n);

  for(int i=0;i<n;i++) dest[i]=(float)(sour[i]);
}

//Change endianess
void doubles_to_floats_changing_endianess(float *dest,double *sour,int n)
{
  char *c;
  char temp;

  if(rank==0 && debug_lvl>1) printf("Converting %d doubles to floats changing endianess\n",n);

  for(int i=0;i<n;i++)
    {
      dest[i]=(float)(sour[i]);

      c=(char*)(dest+i);

      temp=c[3];
      c[3]=c[0];
      c[0]=temp;

      temp=c[2];
      c[2]=c[1];
      c[1]=temp;
    }
}
