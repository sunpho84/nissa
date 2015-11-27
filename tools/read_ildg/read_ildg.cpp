#include "nissa.hpp"

using namespace nissa;

uint64_t max_length=16*1024;

int main(int narg,char **arg)
{
  init_nissa(narg,arg);
  
  //check arguments
  if(narg<2) crash("Use %s file [max_size_to_print=%d]",arg[0],max_length);
  if(narg>=3) max_length=atoi(arg[2]);
  
  //open
  ILDG_File fin=ILDG_File_open_for_read(arg[1]);
  
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
	  
	  //force terminate
	  data[size]='\0';
	  
	  //print and free
	  master_printf("%s\n\n",data);
	  delete[] data;
	}
      else
	{
	  master_printf("skipping the record\n");
	  ILDG_File_skip_record(fin,head);
	}
    }
  
  //close
  ILDG_File_close(fin);
  
  return 0;
}
