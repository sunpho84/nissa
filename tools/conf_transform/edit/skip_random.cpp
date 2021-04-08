#include "nissa.hpp"

using namespace nissa;

int main(int narg,char **arg)
{
  //basic mpi initialization
  init_nissa(narg,arg);
  
  if(narg<7) crash("use: %s T X Y Z file_in file_out",arg[0]);

  for(int mu=0;mu<4;mu++) glbSize[mu]=atoi(arg[mu+1]);

  //Init the MPI grid 
  init_grid(0,0);

  //////////////////////////// read the conf /////////////////////////////

  quad_su3 *conf=nissa_malloc("conf",(locVol+bord_vol).nastyConvert(),quad_su3);  
  //messages
  ILDG_message mess;
  ILDG_message_init_to_last(&mess);
  read_ildg_gauge_conf(conf,arg[5],&mess);
  
  for(ILDG_message *cur_mess=&mess;cur_mess->is_last==false;cur_mess=cur_mess->next)
    if(glb_rnd_gen_inited==0 && strcasecmp(cur_mess->name,"RND_gen_status")==0)
      {
	//grab seed
	start_glb_rnd_gen(cur_mess->data);
	
	//skip
	for(int iskip=0;iskip<10;iskip++) rnd_get_unif(&glb_rnd_gen,0,1);
	
	//store
	free(cur_mess->data);
	char text[1024];
	convert_rnd_gen_to_text(text,&glb_rnd_gen,1024);
	int length=strlen(text)+1;
	cur_mess->data_length=length;
	cur_mess->data=(char*)malloc(length);
	memcpy(cur_mess->data,text,length);
      }
  
  //if message with string not found start from input seed
  if(glb_rnd_gen_inited==0) crash("RND_gen status not found inside conf, starting from input passed seed");
  
  //////////////////////////// write the conf ////////////////////////////
  
  write_ildg_gauge_conf(arg[6],conf,64,&mess);
  nissa_free(conf);
  
  ///////////////////////////////////////////

  close_nissa();

  return 0;
}
