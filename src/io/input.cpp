#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <mpi.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <dirent.h>
#include <filesystem>

#define EXTERN_INPUT
# include "input.hpp"

#include "base/bench.hpp"
#include "base/debug.hpp"
#ifdef USE_TMLQCD
 #include "base/tmLQCD_bridge.hpp"
#endif
#ifdef USE_DDALPHAAMG
 #include "base/DDalphaAMG_bridge.hpp"
#endif
#ifdef USE_QUDA
 #include "base/quda_bridge.hpp"
#endif
#include "base/multiGridParams.hpp"
#include "base/vectors.hpp"
#include "communicate/communicate.hpp"
#include "eigenvalues/eigenvalues.hpp"
#include "geometry/geometry_eo.hpp"
#include "geometry/geometry_lx.hpp"
#include "io/ILDG_File.hpp"
#include "new_types/high_prec.hpp"
#include "new_types/su3.hpp"
#include "routines/ios.hpp"

#include <fcntl.h>
#include <unistd.h>

namespace nissa
{
  //touch a file
  void file_touch(std::string path)
  {
    if(rank==0)
      {
	FILE *f=fopen(path.c_str(),"w");
	if(f!=NULL)
	  {
	    MASTER_PRINTF("File %s created\n",path.c_str());
	    fclose(f);
	  }
	else
	  CRASH("Unable to touch file: %s",path.c_str());
      }
  }
  
  //lock a file
  int file_lock(std::string path)
  {
    int f=0;
    if(rank==0)
      {
	//open the file descriptor
	f=open(path.c_str(),O_RDWR);
	if(f) MASTER_PRINTF("File %s opened\n",path.c_str());
	
	//set a lock and check it
	int status=lockf(f,F_TLOCK,0);
	if(not status)
	  {
	    close(f);
	    f=0;
	    MASTER_PRINTF("Unable to set the lock on %s\n",path.c_str());
	  }
      }
    
    //broadcast
    MPI_Bcast(&f,1,MPI_INT,0,MPI_COMM_WORLD);
    
    return f;
  }
  
  //unlock a file
  int file_unlock(int f)
  {
    if(rank==0)
      {
	int status=lockf(f,F_ULOCK,0);
	if(not status) MASTER_PRINTF("Unable to unset the lock on file descriptor %d\n",f);
	else close(f);
      }
    
    //broadcast
    MPI_Bcast(&f,1,MPI_INT,0,MPI_COMM_WORLD);
    
    return f;
  }
  
  //check if a file exists
  int fileExists(const std::string& path)
  {
    VERBOSITY_LV3_MASTER_PRINTF("Checking if file \"%s\" exists\n",path.c_str());
    
    int status=1;
    
    if(is_master_rank())
      {
	if(FILE *f=fopen(path.c_str(),"r"))
	  {
	    VERBOSITY_LV3_MASTER_PRINTF("File '%s' exists!\n",path.c_str());
	    status=1;
	    fclose(f);
	  }
	else
	  {
	    VERBOSITY_LV3_MASTER_PRINTF("File '%s' do not exist!\n",path.c_str());
	    status=0;
	  }
      }
    
    return broadcast(status);
  }
  
  //return 0 if the dir do not exists, 1 if exists, -1 if exist but is not a directory
  int dir_exists(std::string path)
  {
    int exists;
    
    if(rank==0)
      {
	DIR *d=opendir(path.c_str());
	exists=(d!=NULL);
	
	if(exists)
	  {
	    VERBOSITY_LV2_MASTER_PRINTF("Directory \"%s\" is present\n",path.c_str());
	    closedir(d);
	  }
	else VERBOSITY_LV2_MASTER_PRINTF("Directory \"%s\" is not present\n",path.c_str());
      }
    
    MPI_Bcast(&exists,1,MPI_INT,0,MPI_COMM_WORLD);
    
    return exists;
  }
  
  //open an input file
  void open_input(std::string input_path)
  {
    if(input_global!=NULL) input_global_stack.push_back(input_global);
    
    input_global=open_file(input_path.c_str(),"r");
  }
  
  //close the input file
  void close_input()
  {
    if(rank==0)
      {
	if(input_global_stack.size()==0 and input_global==NULL) CRASH("No input file open");
	fclose(input_global);
      }
    
    if(input_global_stack.size())
      {
	input_global=input_global_stack.back();
	input_global_stack.pop_back();
      }
    else input_global=NULL;
  }
  
  //read a token from file
  int read_next_token(char *tok)
  {
    //avoid warnings
    int ok=1,len=0;
    tok[0]='\0';
    
    if(rank==0)
      {
	if(feof(input_global)) CRASH("reached EOF while scanning input file");
	ok=fscanf(input_global,"%s",tok);
	len=strlen(tok)+1;
      }
    
    MPI_Bcast(&ok,1,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(&len,1,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(tok,len,MPI_BYTE,0,MPI_COMM_WORLD);
    
    return ok;
  }
  
  //check whether a token starts or end a comment
  int check_tok_starts_comment(char *tok)
  {return strncasecmp(tok,"/*",2)==0;}
  int check_tok_ends_comment(char *tok)
  {return strncasecmp(tok+std::max(0lu,strlen(tok)-2),"*/",2)==0;}
  
  //read up to "*/"
  void read_up_to_end_of_comment()
  {
    char tok[1024];
    
    do
      {
	int ok=read_next_token(tok);
	if(ok!=1) CRASH("reached end of file without finding end of comment");
      }
    while(check_tok_ends_comment(tok)==0);
  }
  
  //read variable skipping comments
  int read_var_catcherr(char *out,const char *par,int size_of)
  {
    char tok[1024];
    int ok,comment_found;
    
    do
      {
	//read the token
	ok=read_next_token(tok);
	
	//check if token comment
	if(check_tok_starts_comment(tok))
	  {
	    comment_found=1;
	    //if the comments do not ends itself, read up to finding its end
	    if(!check_tok_ends_comment(tok)) read_up_to_end_of_comment();
	  }
	else comment_found=0;
      }
    while(comment_found);
    
    //parse the token
    if(ok==1)
      ok=(sscanf(tok,par,out)>0);
    
    return ok;
  }
  
  void read_var(char *out,const char *par,int size_of)
  {if(!read_var_catcherr(out,par,size_of)) CRASH("Couldn't read type %s from input file!!!",par);}
  
  //Read an integer from the file
  void read_int(int *out)
  {read_var((char*)out,"%d",sizeof(int));}
  
  //Read a long integer from the file
  void read_int64_t(int64_t *out)
  {read_var((char*)out,"%ld",sizeof(int64_t));}
  
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
    
    if(strcasecmp(exp_str,obt_str)!=0) CRASH("Error, expected '%s' in input file, obtained: '%s'",exp_str,obt_str);
  }
  
  //Read an integer checking the tag
  void read_str_int(const char *exp_str,int *in)
  {
    expect_str(exp_str);
    read_int(in);
    
    VERBOSITY_LV1_MASTER_PRINTF("Read variable '%s' with value: %d\n",exp_str,(*in));
  }
  
  //Read 4 doubles checking the tag
  void read_str_momentum_t(const char *exp_str,
			   Momentum& in)
  {
    expect_str(exp_str);
    for(int i=0;i<4;i++)
      read_double(&in[i]);
    
    VERBOSITY_LV1_MASTER_PRINTF("Read variable '%s' with value: %.16lg %.16lg %.16lg %.16lg\n",exp_str,in[0],in[1],in[2],in[3]);
  }
  
  //Read a double checking the tag
  void read_str_double(const char *exp_str,double *in)
  {
    expect_str(exp_str);
    read_double(in);
    
    VERBOSITY_LV1_MASTER_PRINTF("Read variable '%s' with value: %.16lg\n",exp_str,(*in));
  }
  
  //Read a string checking the tag
  void read_str_str(const char *exp_str,char *in,int length)
  {
    expect_str(exp_str);
    read_str(in,length);
    
    VERBOSITY_LV1_MASTER_PRINTF("Read variable '%s' with value: %s\n",exp_str,in);
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
  
  //Read a list of triples of var, allocating the list
  void read_list_of_var_triples(const char *tag,int *nentries,char **list1,char **list2,char **list3,int size_of_el,const char *par)
  {
    read_str_int(tag,nentries);
    (*list1)=(char*)malloc((*nentries)*size_of_el);
    (*list2)=(char*)malloc((*nentries)*size_of_el);
    (*list3)=(char*)malloc((*nentries)*size_of_el);
    for(int ientr=0;ientr<(*nentries);ientr++)
      {
	char *in1=(*list1)+ientr*size_of_el;
	char *in2=(*list2)+ientr*size_of_el;
	char *in3=(*list3)+ientr*size_of_el;
	read_var(in1,par,size_of_el);
	read_var(in2,par,size_of_el);
	read_var(in3,par,size_of_el);
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
    
    VERBOSITY_LV1_MASTER_PRINTF("List of %s:\t",tag);
    for(int ientr=0;ientr<(*nentries);ientr++) VERBOSITY_LV1_MASTER_PRINTF("%.16lg\t",(*list)[ientr]);
    VERBOSITY_LV1_MASTER_PRINTF("\n");
  }
  void read_list_of_double_pairs(const char *tag,int *nentries,double **list1,double **list2)
  {
    read_list_of_var_pairs(tag,nentries,(char**)list1,(char**)list2,sizeof(double),"%lg");
    
    VERBOSITY_LV1_MASTER_PRINTF("List of pairs of %s:\t",tag);
    for(int ientr=0;ientr<(*nentries);ientr++) VERBOSITY_LV1_MASTER_PRINTF("%.16lg %.16lg\t",(*list1)[ientr],(*list2)[ientr]);
    VERBOSITY_LV1_MASTER_PRINTF("\n");
  }
  void read_list_of_double_triples(const char *tag,int *nentries,double **list1,double **list2,double **list3)
  {
    read_list_of_var_triples(tag,nentries,(char**)list1,(char**)list2,(char**)list3,sizeof(double),"%lg");
    
    VERBOSITY_LV1_MASTER_PRINTF("List of triples of %s:\t",tag);
    for(int ientr=0;ientr<(*nentries);ientr++)
      VERBOSITY_LV1_MASTER_PRINTF("%.16lg %.16lg %.16lg\t",(*list1)[ientr],(*list2)[ientr],(*list3)[ientr]);
    VERBOSITY_LV1_MASTER_PRINTF("\n");
  }
  
  //read a list of int
  void read_list_of_ints(const char *tag,int *nentries,int **list)
  {
    read_list_of_var(tag,nentries,(char**)list,sizeof(int),"%d");
    
    VERBOSITY_LV1_MASTER_PRINTF("List of %s:\t",tag);
    for(int ientr=0;ientr<(*nentries);ientr++) VERBOSITY_LV1_MASTER_PRINTF("%d\t",(*list)[ientr]);
    VERBOSITY_LV1_MASTER_PRINTF("\n");
  }
  
  //read a list of chars
  void read_list_of_chars(const char *tag,int *nentries,char ***list,int nchar_per_entry)
  {
    read_str_int(tag,nentries);
    (*list)=(char**)malloc((*nentries)*sizeof(char*));
    for(int ientr=0;ientr<(*nentries);ientr++)
      {
	(*list)[ientr]=(char*)malloc(nchar_per_entry);
	read_str((*list)[ientr],nchar_per_entry);
      }
    
    VERBOSITY_LV1_MASTER_PRINTF("List of %s:\t",tag);
    for(int ientr=0;ientr<(*nentries);ientr++) VERBOSITY_LV1_MASTER_PRINTF("%s\t",(*list)[ientr]);
    VERBOSITY_LV1_MASTER_PRINTF("\n");
  }
  
  namespace
  {
    /// Ignore unknwon parameters
    int ignoreUnknownParameters;
    
    /// Hold the type tag
    struct triple_tag
    {
      std::string name;
      
      std::string type;
      
      size_t size;
      
      void *pointer;
      
      void read()
      {
	read_var((char*)pointer,type.c_str(),size);
      }
      
      void write()
      {
	if(type=="%d")
	  {
	    VERBOSITY_LV1_MASTER_PRINTF("%d",*((int*)pointer));
	    return;
	  }
	CRASH("unkwnon how to print %s",type.c_str());
      }
      
      const std::string get_tag(int &a)
      {
	return "%d";
      }
      
      template <class T>
      triple_tag(const std::string& name,
		 T &val) :
	name(name),
	type(get_tag(val)),
	size(sizeof(T)),
	pointer(&val)
      {
      }
    };
  }
  
  //read the nissa configuration file
  void read_nissa_config_file()
  {
    char path[1024]="nissa_config";
    
    std::vector<triple_tag> tags;
    tags.push_back(triple_tag("ignore_unknown_parameters",      ignoreUnknownParameters));
    tags.push_back(triple_tag("prepend_time",                   prepend_time));
    tags.push_back(triple_tag("verbosity_lv",                   verbosity_lv));
    tags.push_back(triple_tag("use_128_bit_precision",          use_128_bit_precision));
    tags.push_back(triple_tag("check_inversion_residue",        check_inversion_residue));
    tags.push_back(triple_tag("inversion_residue_heavy_qualify_odg",inversion_residue_heavy_qualify_odg));
    tags.push_back(triple_tag("inversion_residue_threshold_odg",inversion_residue_threshold_odg));
    tags.push_back(triple_tag("use_eo_geom",                    use_eo_geom));
    tags.push_back(triple_tag("use_async_communications",       use_async_communications));
    tags.push_back(triple_tag("warn_if_not_disallocated",       warn_if_not_disallocated));
    tags.push_back(triple_tag("set_t_nranks",                   fix_nranks[0]));
    tags.push_back(triple_tag("set_x_nranks",                   fix_nranks[1]));
    tags.push_back(triple_tag("set_y_nranks",                   fix_nranks[2]));
    tags.push_back(triple_tag("set_z_nranks",                   fix_nranks[3]));
    tags.push_back(triple_tag("ignore_ILDG_magic_number",       ignore_ILDG_magic_number));
    tags.push_back(triple_tag("fast_read_write_vectors",        fast_read_write_vectors));
    tags.push_back(triple_tag("perform_benchmark",              perform_benchmark));
#if HIGH_PREC_TYPE==GMP_HIGH_PREC
    tags.push_back(triple_tag("mpf_precision",                  mpf_precision));
#endif
#ifdef USE_TMLQCD
    tags.push_back(triple_tag("use_tmLQCD",                     use_tmLQCD));
#endif
#if defined USE_DDALPHAAMG or defined USE_QUDA
    tags.push_back(triple_tag("use_DDalphaAMG",                 multiGrid::use_multiGrid));
    tags.push_back(triple_tag("use_multigrid",                  multiGrid::use_multiGrid));
    tags.push_back(triple_tag("use_deflated_solver",            multiGrid::use_deflated_solver));
#endif
#ifdef USE_QUDA
    tags.push_back(triple_tag("use_QUDA",                       use_quda));
#endif
#ifdef USE_PARPACK
    tags.push_back(triple_tag("use_parpack",                    use_parpack));
#endif
    
    if(fileExists(path))
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
		size_t itag=0;
		while(itag<tags.size() && tag!=tags[itag].name) itag++;
		
		//check if tag found
		if(itag==tags.size())
		  if(ignoreUnknownParameters)
		    {
		      WARNING("unkwnown parameter '%s'\n",tag);
		      char tmp[100];
		      read_var_catcherr(tmp,"%s",100);
		      WARNING(" skipped value '%s' for unknown parameter \'%s\'\n",tmp,tag);
		    }
		  else
		    {
		      MASTER_PRINTF("Available parameters:\n");
		      for(size_t itag=0;itag<tags.size();itag++)
			MASTER_PRINTF(" %s\n",tags[itag].name.c_str());
		      
		      CRASH("unkwnown parameter '%s', set \"ignore_unknown_parameters to 1\" to bypass",tag);
		    }
		else
		  {
		    //read the tag
		    tags[itag].read();
		    VERBOSITY_LV1_MASTER_PRINTF("Read parameter '%s' with value ",tag);
		    tags[itag].write();
		    VERBOSITY_LV1_MASTER_PRINTF("\n");
		  }
	      }
	    else
	      MASTER_PRINTF("Finished reading the file '%s'\n",path);
	  }
	while(nr==1);
	
	close_input();
      }
    else MASTER_PRINTF("No 'nissa_config' file present, using standard configuration\n");
    
    VERBOSITY_LV1_MASTER_PRINTF("Configuration:\n");
    for(size_t itag=0;itag<tags.size();itag++)
      {
	VERBOSITY_LV1_MASTER_PRINTF(" %s=",tags[itag].name.c_str());
	tags[itag].write();
	VERBOSITY_LV1_MASTER_PRINTF("\n");
      }
  }
}
