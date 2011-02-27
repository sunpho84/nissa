#pragma once
#include <stdio.h>

void expect_string_from_file(FILE *fin,const char *exp)
{
  char read[1024];
  int nscan=fscanf(fin,"%s",read);
  if(nscan!=1||strcmp(exp,read)!=0)
    {
      if(nscan==EOF) fprintf(stderr,"Error,reached file end while waiting for '%s'\n",exp);
      else           fprintf(stderr,"Error, read '%s' instead than '%s'\n",read,exp);
      exit(1);
    }
}

void read_formatted_from_file(char *out,FILE *fin,const char *what,const char *varname)
{
  int nscan=fscanf(fin,what,out);
  if(nscan!=1)
    {
      if(nscan==EOF) fprintf(stderr,"Error,reached file end while reading '%s'\n",varname);
      else           fprintf(stderr,"Error, not enough data while reading '%s'\n",varname);
      exit(1);
    }
}

void read_formatted_from_file_expecting(char *out,FILE *fin,const char *what,const char *varname)
{
  expect_string_from_file(fin,varname);
  read_formatted_from_file(out,fin,what,varname);
}

FILE *open_file(const char* path,const char* mod)
{
  FILE *out=fopen(path,mod);

  if(out==NULL)
    {
      fprintf(stderr,"Error opening file '%s' in mode '%s'\n",path,mod);
      exit(1);
    }
  
  return out;
}

