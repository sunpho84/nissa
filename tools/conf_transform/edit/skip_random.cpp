#include "nissa.hpp"

using namespace nissa;

int main(int narg,char **arg)
{
  //basic mpi initialization
  init_nissa(narg,arg);
  
  if(narg<5) crash("use: %s L T file_in file_out",arg[0]);

  int L=atoi(arg[1]);
  int T=atoi(arg[2]);

  //Init the MPI grid 
  init_grid(T,L);

  //////////////////////////// read the conf /////////////////////////////

  quad_su3 *conf=nissa_malloc("conf",loc_vol+bord_vol,quad_su3);  
  //messages
  ILDG_message mess;
  ILDG_message_init_to_last(&mess);
  read_ildg_gauge_conf(conf,arg[3],&mess);
  
  for(ILDG_message *cur_mess=&mess;cur_mess->is_last==false;cur_mess=cur_mess->next)
    if(glb_rnd_gen_inited==0 && strcasecmp(cur_mess->name,"RND_gen_status")==0)
      {
	//grab seed
	start_loc_rnd_gen(cur_mess->data);
	
	//skip
	for(int iskip=0;iskip<10;iskip++) rnd_get_unif(&glb_rnd_gen,0,1);
	
	//store
	convert_rnd_gen_to_text(cur_mess->data,&glb_rnd_gen);
      }
  
  //if message with string not found start from input seed
  if(glb_rnd_gen_inited==0) crash("RND_gen status not found inside conf, starting from input passed seed");
  
  //////////////////////////// write the conf ////////////////////////////
  
  write_ildg_gauge_conf(arg[4],conf,64,&mess);
  nissa_free(conf);
  
  ///////////////////////////////////////////

  close_nissa();

  return 0;
}
