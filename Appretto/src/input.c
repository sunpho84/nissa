#pragma once

#include <strings.h>

FILE *input_global;

void open_input(char *input_path)
{
  if(rank==0)
    {
      input_global=fopen(input_path,"r");
      if(input_global==NULL)
	{
	  fprintf(stderr,"File '%s' not found\n",input_path);
	  fflush(stderr);
	  MPI_Abort(MPI_COMM_WORLD,1);
	}

      if(debug) printf("File '%s' opened\n",input_path);
    }	
}

void close_input()
{
  if(rank==0) fclose(input_global);
}

//Read an integer from the file
void read_int(int *in)
{
  if(rank==0)
    {
      int ok=fscanf(input_global,"%d",in);
      if(!ok)
	{
	  fprintf(stderr,"Couldn't read from input file!!!\n");
	  fflush(stderr);
	  MPI_Abort(MPI_COMM_WORLD,1);
	}
    }
  MPI_Bcast(in,1,MPI_INT,0,MPI_COMM_WORLD);
}

//Read a double from the file
void read_double(double *in)
{
  if(rank==0)
    {
      int ok=fscanf(input_global,"%lg",in);
      if(!ok)
	{
	  fprintf(stderr,"Couldn't read from input file!!!\n");
	  fflush(stderr);
	  MPI_Abort(MPI_COMM_WORLD,1);
	}
    }
  MPI_Bcast(in,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
}

//Read a string from the file
void read_str(char *str,int length)
{
  if(rank==0)
    {
      int ok=fscanf(input_global,"%s",str);
      if(!ok)
	{
	  fprintf(stderr,"Couldn't read from input file!!!\n");
	  fflush(stderr);
	  MPI_Abort(MPI_COMM_WORLD,1);
	}
    }
  MPI_Bcast(str,length,MPI_BYTE,0,MPI_COMM_WORLD);
}

//Read a string from the file and check against the argument
void expect_str(const char *exp_str)
{
  char obt_str[1024];

  read_str(obt_str,1024);
  
  if(strcasecmp(exp_str,obt_str)!=0 && rank==0)
    {
      fprintf(stderr,"Error, expexcted '%s' in input file, obtained: '%s'\n",exp_str,obt_str);
      fflush(stderr);
      MPI_Abort(MPI_COMM_WORLD,1);
    }
}

//Read an integer checking the tag
void read_str_int(const char *exp_str,int *in)
{
  expect_str(exp_str);
  read_int(in);

  if(rank==0 && debug) printf("Read variable '%s' with value: %d\n",exp_str,(*in));
}

//Read a double checking the tag
void read_str_double(const char *exp_str,double *in)
{
  expect_str(exp_str);
  read_double(in);

  if(rank==0 && debug) printf("Read variable '%s' with value: %f\n",exp_str,(*in));
}

//Read a string checking the tag
void read_str_str(const char *exp_str,char *in,int length)
{
  expect_str(exp_str);
  read_str(in,length);

  if(rank==0 && debug) printf("Read variable '%s' with value: %s\n",exp_str,in);
}
