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
	  verbosity_lv2_master_printf("File '%s' exists!\n",path);
          status=1;
          fclose(f);
        }
      else
        {
	  verbosity_lv2_master_printf("File '%s' do not exist!\n",path);
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

//open an input file
void open_input(char *input_path)
{
  if(rank==0)
    {
      input_global=fopen(input_path,"r");
      if(input_global==NULL) crash("File '%s' not found",input_path);
      
      verbosity_lv1_master_printf("File '%s' opened\n",input_path);
    }	
}

//close the input file
void close_input()
{if(rank==0) fclose(input_global);}

//Read a var from the file
int read_var_catcherr(char *out,const char *par,int size_of)
{
  int ok=(rank==0) ? fscanf(input_global,par,out) : 1;
  MPI_Bcast(&ok,1,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(out,size_of,MPI_BYTE,0,MPI_COMM_WORLD);
  return ok;
}

void read_var(char *out,const char *par,int size_of)
{if(!read_var_catcherr(out,par,size_of)) crash("Couldn't read from input file!!!");}

//Read an integer from the file
void read_int(int *out)
{read_var((char*)out,"%d",sizeof(int));}

//Read a double from the file
void read_double(double *out)
{read_var((char*)out,"%lg",sizeof(double));}

//Read a string from the file
void read_str(char *str,int length)
{read_var(str,"%s",length);}

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

//Read a list of var, allocating the list
void read_list_of_var(const char *tag,int *nentries,char **list,int size_of_el,const char *par)
{
  read_str_int(tag,nentries);
  (*list)=(char*)malloc((*nentries)*size_of_el);
  for(int ientr=0;ientr<(*nentries);ientr++)
    {
      char *in=(*list)+ientr*size_of_el;
      read_var(in,par,size_of_el);
    }
}

//Read a list of pairs of var, allocating the list
void read_list_of_var_pairs(const char *tag,int *nentries,char **list1,char **list2,int size_of_el,const char *par)
{
  read_str_int(tag,nentries);
  (*list1)=(char*)malloc((*nentries)*size_of_el);
  (*list2)=(char*)malloc((*nentries)*size_of_el);
  for(int ientr=0;ientr<(*nentries);ientr++)
    {
      char *in1=(*list1)+ientr*size_of_el;
      char *in2=(*list2)+ientr*size_of_el;
      read_var(in1,par,size_of_el);
      read_var(in2,par,size_of_el);
    }
}

//read a list of doubles
void read_list_of_doubles(const char *tag,int *nentries,double **list)
{
  read_list_of_var(tag,nentries,(char**)list,sizeof(double),"%lg");
  
  verbosity_lv1_master_printf("List of %s:\t",tag);
  for(int ientr=0;ientr<(*nentries);ientr++) verbosity_lv1_master_printf("%lg\t",(*list)[ientr]);
  verbosity_lv1_master_printf("\n");
}

//read a list of doubles
void read_list_of_double_pairs(const char *tag,int *nentries,double **list1,double **list2)
{
  read_list_of_var_pairs(tag,nentries,(char**)list1,(char**)list2,sizeof(double),"%lg");
  
  verbosity_lv1_master_printf("List of pairs of %s:\t",tag);
  for(int ientr=0;ientr<(*nentries);ientr++) verbosity_lv1_master_printf("%lg %lg\t",(*list1)[ientr],(*list2)[ientr]);
  verbosity_lv1_master_printf("\n");
}

//read a list of int
void read_list_of_ints(const char *tag,int *nentries,int **list)
{
  read_list_of_var(tag,nentries,(char**)list,sizeof(int),"%d");
  
  verbosity_lv1_master_printf("List of %s:\t",tag);
  for(int ientr=0;ientr<(*nentries);ientr++) verbosity_lv1_master_printf("%d\t",(*list)[ientr]);
  verbosity_lv1_master_printf("\n");
}

//read the nissa configuration file
void read_nissa_config_file()
{
  char path[1024]="nissa_config";
  
  const int navail_tag=2;
  char tag_name[2][100]={"verbosity_lv","use_128_bit_precision"};
  char *tag_addr[2]={(char*)&nissa_verbosity,(char*)&nissa_use_128_bit_precision};
  char tag_type[2][3]={"%d","%d"};
  char tag_size[2]={4,4};
  
  if(file_exists(path))
    {
      open_input(path);
      int nr;
      
      do
        {
          char tag[100];
	  nr=read_var_catcherr(tag,"%s",100);
	  
	  if(nr>=1)
	    {
	      //find the tag
	      int itag=0;
	      while(itag<navail_tag && strcasecmp(tag,tag_name[itag])!=0) itag++;
	      
	      //check if tag found
	      if(itag==navail_tag) crash("unkwnown parameter '%s'",tag);
	      
	      //read the tag
	      read_var(tag_addr[itag],tag_type[itag],tag_size[itag]);
	      verbosity_lv1_master_printf("Read parameter '%s' with value ",tag);
	      verbosity_lv1_master_printf(tag_type[itag],*((int*)tag_addr[itag]));
	      verbosity_lv1_master_printf("\n");
	    }
	  else master_printf("Finished reading the file '%s'\n",path);
	}
      while(nr==1);
      
      close_input();
    }
  else master_printf("No 'nissa_config' file present, using standard configuration\n");
  
  verbosity_lv1_master_printf("Configuration:\n");
  verbosity_lv1_master_printf(" verbosity_lv=%d\n",nissa_verbosity);
  verbosity_lv1_master_printf(" use_128_bit_precision=%d\n",nissa_use_128_bit_precision);
}
