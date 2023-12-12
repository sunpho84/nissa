#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#define EXTERN_IOS
# include "ios.hpp"

#include <math.h>
#include <mpi.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <string.h>
#include <sys/stat.h>

#include <base/debug.hpp>
#include <routines/mpi_routines.hpp>
#include <threads/threads.hpp>

namespace nissa
{
  //return the number of occurency of "sub" inside "str"
  int count_substrings(const char *str,const char *sub)
  {
    int rc=0;
    for(str=(str!=NULL)?strstr(str,sub):str;str!=NULL;str=strstr(str+1,sub)) rc++;
    
    return rc;
  }
  
  //get a unique temporary filename
  void master_get_temp_file(FILE *&fout,std::string &prefix)
  {
    //prepare the full name
    char *buffer=strdup((prefix+"XXXXXX").c_str());
    
    //open on the master
    int fd=-1;
    if(isMasterRank())
      {
	fd=mkstemp(buffer);
	if(fd==-1) crash("failed to open a temporary file with prefix %s",prefix.c_str());
	fout=fdopen(fd,"w");
      }
    
    //broadcast name and copy in prefix
    crash("");
    // mpiBcast(buffer);
    // MPI_Bcast(buffer,strlen(buffer),MPI_CHAR,rank,MPI_COMM_WORLD);
    prefix=buffer;
    
    free(buffer);
  }
  
  //only master rank and thread print
  int master_fprintf(FILE *stream,const char *format,...)
  {
    int ret=0;
    
    static bool print_time=true;
    
    if(prepend_time and print_time and isMasterRank())
      ret+=fprintf(stream,"%lg s:\t",take_time());
    
    va_list ap;
    va_start(ap,format);
    if(isMasterRank())
      ret=vfprintf(stream,format,ap);
    va_end(ap);
    
    print_time=(format[strlen(format)-1]=='\n');
    
    return ret;
  }
  
  //print a number in some kind of familiar units
  void fprintf_friendly_units(FILE *fout,uint64_t quant,uint64_t orders,const char *units)
  {
    const char units_prefix[6][2]={"","K","M","G","T","P"};
    int iord=0;
    double temp=quant;
    
    while(temp>orders)
      {
	temp/=orders;
	iord++;
      }
    
    master_fprintf(fout,"%d %s%s",(int)(temp+0.5),units_prefix[iord],units);
  }
  void fprintf_friendly_filesize(FILE *fout,uint64_t quant)
  {fprintf_friendly_units(fout,quant,1024,"Bytes");}
  
  //create a dir
  int create_dir(std::string path)
  {
    umask(0);
    const int res=getExecOnMasterRank(mkdir,path.c_str(),0775);
    if(res!=0)
      master_printf("Warning, failed to create dir %s, returned %d. Check that you have permissions and that parent dir exists.\n",path.c_str(),res);
    
    return res;
  }
  
  //copy a file
  int cp(std::string path_out,std::string path_in)
  {
    int rc=0;
    if(isMasterRank())
      {
	char command[1024];
	sprintf(command,"cp %s %s",path_in.c_str(),path_out.c_str());
	rc=system(command);
	if(rc!=0) crash("cp failed!");
      }
    
    return broadcast(rc);
  }
  
  //pass to the folder
  int cd(std::string path)
  {
    int rc=0;
    if(isMasterRank())
      {
	char command[1024];
	sprintf(command,"cd %s",path.c_str());
	rc=system(command);
	if(rc!=0) crash("cd to %s failed!",path.c_str());
      }
    
    return broadcast(rc);
  }
  
  //Open a file checking it
  FILE *open_file(std::string outfile,const char *mode,const MpiRank& ext_rank)
  {
    FILE *fout=NULL;
    
    if(ext_rank==-1 or thisRank==ext_rank)
      {
	if(outfile=="-") fout=stdout;
	else
	  {
	    fout=fopen(outfile.c_str(),mode);
	    if(fout==NULL) crash("Couldn't open file: \"%s\" with mode: \"%s\"",outfile.c_str(),mode);
	  }
      }
    
    return fout;
  }
  
  //Open a text file for output
  FILE* open_text_file_for_output(std::string outfile)
  {return open_file(outfile,"w");}
  
  //Open a text file for input
  FILE* open_text_file_for_input(std::string infile)
  {return open_file(infile,"r");}
  
  //close an open file
  void close_file(FILE *file)
  {if(isMasterRank() && file!=stdout) fclose(file);}
  
  //check if a file exists
  int file_exists(std::string path)
  {
    int status=1;
    
    if(isMasterRank()==0)
      {
	FILE *f=fopen(path.c_str(),"r");
	if(f!=NULL)
	  {
	    verbosity_lv3_master_printf("File '%s' exists!\n",path.c_str());
	    status=1;
	    fclose(f);
	  }
	else
	  {
	    verbosity_lv3_master_printf("File '%s' do not exist!\n",path.c_str());
	    status=0;
	  }
      }
    
    //broadcast the result
    MPI_Bcast(&status,1,MPI_INT,0,MPI_COMM_WORLD);
    
    return status;
  }
  
  //count the number of lines in a file
  int count_file_lines(std::string path)
  {
    //return -1 if file does not exist
    if(!file_exists(path)) return -1;
    
    //scan the file
    FILE *fin=open_text_file_for_input(path);
    int n=0;
    if(isMasterRank())
      {
	int ch;
	while(EOF!=(ch=getchar())) if (ch=='\n') n++;
      }
    
    //close file and broadcast n
    close_file(fin);
    
    return broadcast(n);
  }
  
  //get the size of a file
  int get_file_size(std::string path)
  {
    //return -1 if file does not exist
    if(!file_exists(path)) return -1;
    
    //scan the file
    FILE *fin=open_text_file_for_input(path);
    int file_size=0;
    if(isMasterRank())
      {
	if(fseek(fin,0,SEEK_END)) crash("while seeking");
	file_size=ftell(fin);
      }
    
    return broadcast(file_size);
  }
  
  //take the last characters of the passed string
  void take_last_characters(char *out,const char *in,int size)
  {
    int len=strlen(in)+1;
    int copy_len=(len<=size)?len:size;
    const char *str_init=(len<=size)?in:in+len-size;
    ::memcpy(out,str_init,copy_len);
  }
  
  //combine arguments in a single string
  std::string combine(const char *format,...)
  {
    char buffer[1024];
    va_list args;
    
    va_start(args,format);
    vsnprintf(buffer,1024,format,args);
    va_end(args);
    
    return std::string(buffer);
  }
}
