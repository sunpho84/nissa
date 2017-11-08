#include "nissa.hpp"

using namespace nissa;

uint64_t max_length=16*1024;
char *file_pattern=NULL;

int main(int narg,char **arg)
{
  init_nissa(narg,arg);
  
  //check arguments
  if(narg<2) crash("Use %s file [max_size_to_print=%d] [file_pattern]",arg[0],max_length);
  if(narg>=3) max_length=atoi(arg[2]);
  if(narg>=4) file_pattern=arg[3];
  
  //open
  ILDG_File fin=ILDG_File_open_for_read(arg[1]);

  int irec=0;
  while(!ILDG_File_reached_EOF(fin))
    {
      //get header
      ILDG_header head=ILDG_File_get_next_record_header(fin);
      
      //read size
      uint64_t size=head.data_length;
      master_printf("Found record: %s, length: %lu\n",head.type,size);
      
      //read if not too big
      if(size<max_length)
	{
	  //allocate and read
	  char *data=new char[size+1];
	  ILDG_File_read_all(data,fin,size);
	  
	  if(file_pattern==NULL)
	    {
	      //force terminate
	      data[size]='\0';
	      
	      //print
	      master_printf("%s\n\n",data);
	    }
	  else
	    {
	      FILE *fout=open_file(combine("%s_rec_%d_%s",file_pattern,irec,head.type),"w");
	      if(fwrite(data, sizeof(char),size,fout)!=size) crash("writing record %s",head.type);
	      close_file(fout);
	    }
	  delete[] data;
	}
      else
	{
	  master_printf("skipping the record\n");
	  ILDG_File_skip_record(fin,head);
	}
      
      irec++;
    }
  
  //close
  ILDG_File_close(fin);
  
  return 0;
}
