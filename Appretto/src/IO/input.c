#pragma once

FILE *input_global;

int file_exist(const char *path)
{
  int status=1;
  
  if(rank==0)
    {
      FILE *f=fopen(path,"r");
      if(f!=NULL)
        {
          status=1;
          fclose(f);
        }
      else status=0;
    }
  
  MPI_Bcast(&status,1,MPI_INT,0,MPI_COMM_WORLD);
  return status;
}

void open_input(char *input_path)
{
  if(rank==0)
    {
      input_global=fopen(input_path,"r");
      if(input_global==NULL) crash("File '%s' not found",input_path);

      if(debug) printf("File '%s' opened\n",input_path);
    }	
}

void close_input()
{
  if(rank==0) fclose(input_global);
}

//Read a var from the file
void read_var(char *in,const char *par,int size_of)
{
  int ok=(rank==0) ? fscanf(input_global,par,in) : 1;
  if(!ok) crash("Couldn't read from input file!!!");
  
  MPI_Bcast(in,size_of,MPI_BYTE,0,MPI_COMM_WORLD);
}

//Read an integer from the file
void read_int(int *in)
{read_var((char*)in,"%d",sizeof(int));}

//Read a double from the file
void read_double(double *in)
{read_var((char*)in,"%lg",sizeof(double));}

//Read a string from the file
void read_str(char *str,int length)
{read_var((char*)str,"%s",length);}

//Read a string from the file and check against the argument
void expect_str(const char *exp_str)
{
  char obt_str[1024];

  read_str(obt_str,1024);
  
  if(strcasecmp(exp_str,obt_str)!=0) crash("Error, expexcted '%s' in input file, obtained: '%s'",exp_str,obt_str);
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

  if(rank==0 && debug) printf("Read variable '%s' with value: %g\n",exp_str,(*in));
}

//Read a string checking the tag
void read_str_str(const char *exp_str,char *in,int length)
{
  expect_str(exp_str);
  read_str(in,length);

  if(rank==0 && debug) printf("Read variable '%s' with value: %s\n",exp_str,in);
}

//Read a list of double and its length, allocate the list
void read_list_of_var(char *tag,int *nentries,char **list,int size_of_el,const char *par)
{
  read_str_int(tag,nentries);
  (*list)=(char*)malloc((*nentries)*size_of_el);
  for(int ientr=0;ientr<(*nentries);ientr++)
    {
      char *in=(*list)+ientr*size_of_el;
      read_var(in,par,size_of_el);
    }
}

//read a list of doubles
void read_list_of_doubles(char *tag,int *nentries,double **list)
{
  read_list_of_var(tag,nentries,(char**)list,sizeof(double),"%lg");
  if(rank==0)
    {
      printf("List of %s:\t",tag);
      for(int ientr=0;ientr<(*nentries);ientr++) printf("%lg\t",(*list)[ientr]);
      printf("\n");
    }
}

//read a list of int
void read_list_of_ints(char *tag,int *nentries,int **list)
{
  read_list_of_var(tag,nentries,(char**)list,sizeof(int),"%d");
  if(rank==0)
    {
      printf("List of %s:\t",tag);
      for(int ientr=0;ientr<(*nentries);ientr++) printf("%d\t",(*list)[ientr]);
      printf("\n");
    }
}
