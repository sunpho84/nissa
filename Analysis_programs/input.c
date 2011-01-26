#pragma once

#include <stdio.h>
#include <stdlib.h>

FILE *input_global;

//Read an integer from the file
void read_int(int *in)
{
  int ok=fscanf(input_global,"%d",in);
  if(!ok)
    {
      fprintf(stderr,"Couldn't read from input file!!!\n");
      exit(1);
    }
}

//Read a double from the file
void read_double(double *in)
{
  int ok=fscanf(input_global,"%lg",in);
  if(!ok)
    {
      fprintf(stderr,"Couldn't read from input file!!!\n");
      exit(1);
    }
}

//Read a string from the file
void read_str(char *str,int length)
{
  int ok=fscanf(input_global,"%s",str);
  if(!ok)
    {
      fprintf(stderr,"Couldn't read from input file!!!\n");
      exit(1);
    }
}

//Read a string from the file and check against the argument
void expect_str(const char *exp_str)
{
  char obt_str[1024];

  read_str(obt_str,1024);
  
  if(strcasecmp(exp_str,obt_str)!=0)
    {
      fprintf(stderr,"Error, expexcted '%s' in input file, obtained: '%s'\n",exp_str,obt_str);
      exit(1);
    }
}

//Read an integer checking the tag
void read_str_int(const char *exp_str,int *in)
{
  expect_str(exp_str);
  read_int(in);
}

//Read a double checking the tag
void read_str_double(const char *exp_str,double *in)
{
  expect_str(exp_str);
  read_double(in);
}

//Read a string checking the tag
void read_str_str(const char *exp_str,char *in,int length)
{
  expect_str(exp_str);
  read_str(in,length);
}

void open_input(char *input_path)
{
  input_global=fopen(input_path,"r");
  if(input_global==NULL)
    {
      fprintf(stderr,"File '%s' not found\n",input_path);
      exit(1);
    }
}
