#include "nissa.hpp"

using namespace nissa;

int main(int narg,char **arg)
{
  init_nissa(narg,arg);
  
  if(nranks>1) crash("cannot run in parallel");
  
  if(narg<7) crash("use: %s L T file_in file_out parity nspinors",arg[0]);
  
  int L=atoi(arg[1]);
  int T=atoi(arg[2]);
  char *pathin=arg[3];
  char *pathout=arg[4];
  int parity=atoi(arg[5]);
  int nspinors=atoi(arg[6]);
  
  //Init the MPI grid
  init_grid(T,L);
  
  ///////////////////////////////////////////
  
  spincolor *sp=nissa_malloc("sp",locVol.nastyConvert(),spincolor);
  int ispinor=0;
  ILDG_File fin=ILDG_File_open_for_read(pathin);
  ILDG_File fout=ILDG_File_open_for_write(pathout);
  do
    {
      master_printf("Searching for spinor: %d/%d\n",ispinor,nspinors);
      
      ILDG_header head=ILDG_File_get_next_record_header(fin);
      master_printf("%s %lld\n",head.type,head.data_length);
      
      if(strcasecmp(head.type,"scidac-binary-data")==0)
	{
	  read_real_vector(sp,fin,head);
	  NISSA_LOC_VOL_LOOP(ivol)
	    if(loclx_parity[ivol.nastyConvert()]!=parity)
	      spincolor_put_to_zero(sp[ivol.nastyConvert()]);
	  write_real_vector(fout,sp,64,head.type);
	  ispinor++;
	}
      else
	{
	  char *mess=(char*)malloc(head.data_length+1);
	  ILDG_File_read_all(mess,fin,head.data_length);
	  ILDG_File_write_record(fout, head.type,mess,head.data_length);
	  
	  //write
	  mess[head.data_length]='\0';
	  master_printf("%s\n================================================\n",mess);
	  free(mess);
	}
    }
  while(ispinor<nspinors);
  
  ILDG_File_close(fin);
  ILDG_File_close(fout);
  
  nissa_free(sp);
  
  return 0;
}
