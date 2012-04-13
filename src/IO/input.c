#pragma once

FILE *input_global;

void file_touch(const char *path)
{
  if(rank==0)
    {
      FILE *f=fopen(path,"w");
      if(f!=NULL)
        {
	  master_printf("File %s created\n",path);
          fclose(f);
	}
      else
	crash("Unable to touch file: %s\n",path);
    }
}

int file_exists(const char *path)
{
  int status=1;
  
  if(rank==0)
    {
      FILE *f=fopen(path,"r");
      if(f!=NULL)
        {
	  master_printf("File %s exists!\n",path);
          status=1;
          fclose(f);
        }
      else
        {
	  master_printf("File %s do not exist!\n",path);
	  status=0;
	}
    }
  
  MPI_Bcast(&status,1,MPI_INT,0,MPI_COMM_WORLD);
  return status;
}

//return 0 if the dir do not exists, 1 if exists, -1 if exist but is not a directory
int dir_exists(char *path)
{
  struct stat st;
  int exists;
  
  if(rank==0)
    {
      exists=(stat(path,&st)==0);

      if(exists) verbosity_lv2_master_printf("Directory \"%s\" is present\n",path);
      else verbosity_lv2_master_printf("Directory \"%s\" is not present\n",path);
    }
  
  MPI_Bcast(&exists,1,MPI_INT,0,MPI_COMM_WORLD);
  
  return exists;
}

void open_input(char *input_path)
{
  if(rank==0)
    {
      input_global=fopen(input_path,"r");
      if(input_global==NULL) crash("File '%s' not found",input_path);
      
      verbosity_lv1_master_printf("File '%s' opened\n",input_path);
    }	
}

void close_input()
{
  if(rank==0) fclose(input_global);
}

//Read a var from the file
int read_var_catcherr(char *in,const char *par,int size_of)
{
  int ok=(rank==0) ? fscanf(input_global,par,in) : 1;
  MPI_Bcast(in,size_of,MPI_BYTE,0,MPI_COMM_WORLD);
  return ok;
}

void read_var(char *in,const char *par,int size_of)
{if(!read_var_catcherr(in,par,size_of)) crash("Couldn't read from input file!!!");}

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

  verbosity_lv1_master_printf("Read variable '%s' with value: %d\n",exp_str,(*in));
}

//Read a double checking the tag
void read_str_double(const char *exp_str,double *in)
{
  expect_str(exp_str);
  read_double(in);

  verbosity_lv1_master_printf("Read variable '%s' with value: %g\n",exp_str,(*in));
}

//Read a string checking the tag
void read_str_str(const char *exp_str,char *in,int length)
{
  expect_str(exp_str);
  read_str(in,length);

  verbosity_lv1_master_printf("Read variable '%s' with value: %s\n",exp_str,in);
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
  
  master_printf("List of %s:\t",tag);
  for(int ientr=0;ientr<(*nentries);ientr++) master_printf("%lg\t",(*list)[ientr]);
  master_printf("\n");
}

//read a list of int
void read_list_of_ints(char *tag,int *nentries,int **list)
{
  read_list_of_var(tag,nentries,(char**)list,sizeof(int),"%d");
  
  master_printf("List of %s:\t",tag);
  for(int ientr=0;ientr<(*nentries);ientr++) master_printf("%d\t",(*list)[ientr]);
  master_printf("\n");
}

//read the nissa configuration file
void read_nissa_config_file()
{
  const int navail_tag=1;
  char tag_name[1][100]={"nissa_verbosity_lv"};
  void *tag_addr[1]={&nissa_verbosity};
  char tag_type[1][3]={"%d"};
  char tag_size[1]={4};
  
  if(file_exists("nissa_config"))
    {
      open_input("nissa_config");
      int ok;
      
      do
        {
          char tag[100];
	  ok=read_str_catherr(tag,100);
	  
	  if(ok)
	    {
	      //find the tag
	      int itag=0;
	      while(itag<nvail_tag && strcasecmp(tag,tag_name[itag])!=0) itag++;
	      
	      //check if tag found
	      if(itag==nvail_tag) crash("unkwnown tag '%s'",tag);
	      
	      //read the tag
	      read_var((char*)tag_addr[itag],tag_type[itag],tag_size[itag]);
	    }
	  else master_printf("Finished reading the file");
	}
      while(ok);
      
      close_input();
    }
  else master_printf("No 'nissa_config' file present, using standard configuration\n");
}
