#include "nissa.hpp"

using namespace nissa;

bool rnd_gen_diff(rnd_gen &a,rnd_gen &b,int tot)
{
  bool diff=(a.idum!=b.idum);
  diff|=(a.idum2!=b.idum2);
  for(int i=0;i<RAN2_NTAB;i++) diff|=(a.iv[i]!=b.iv[i]);
  diff|=(a.iy!=b.iy);
  char tempa[1024],tempb[1024];
  convert_rnd_gen_to_text(tempa,&a,1024);
  convert_rnd_gen_to_text(tempb,&b,1024);
  MASTER_PRINTF("iter %d %s\n",tot,tempa);
  MASTER_PRINTF("iter %d %s\n",tot,tempb);
  MASTER_PRINTF("\n");
  return diff;
}

int main(int narg,char **arg)
{
  //basic mpi initialization
  initNissa(narg,arg);
  
  if(narg<8) CRASH("use: %s T X Y Z file_in file_out ori_seed",arg[0]);
  
  for(int mu=0;mu<4;mu++) glbSize[mu]=atoi(arg[mu+1]);
  
  //Init the MPI grid 
  initGrid(0,0);
  
  //////////////////////////// read the conf /////////////////////////////
  
  quad_su3 *conf=nissa_malloc("conf",locVol+bord_vol,quad_su3);
  //messages
  ILDG_message mess;
  ILDG_message_init_to_last(&mess);
  read_ildg_gauge_conf(conf,arg[5],&mess);
  
  for(ILDG_message *cur_mess=&mess;cur_mess->is_last==false;cur_mess=cur_mess->next)
    if(glb_rnd_gen_inited==0 && strcasecmp(cur_mess->name,"RND_gen_status")==0)
      {
	//grab ref status
	rnd_gen ref;
	convert_text_to_rnd_gen(&ref,cur_mess->data);
	
	//skip
	const int nrew=10;
	rnd_gen list[nrew+1];
	start_rnd_gen(&(list[0]),atoi(arg[7]));
	int cur=0;
	int tot=0;

	while(rnd_gen_diff(list[cur],ref,tot))
	  {
	    int next=(cur+1)%(nrew+1);
	    memcpy(list+next,list+cur,sizeof(rnd_gen));
	    rnd_get_unif(list+next,0,1);
	    cur=next;
	    tot++;
	  }
	MASTER_PRINTF("followed %d extraction\n",tot);
	
	//store
	free(cur_mess->data);
	char text[1024];
	convert_rnd_gen_to_text(text,list+(cur+1)%(nrew+1),1024);
	int length=strlen(text)+1;
	cur_mess->data_length=length;
	cur_mess->data=(char*)malloc(length);
	memcpy(cur_mess->data,text,length);
      }
  
  //////////////////////////// write the conf ////////////////////////////
  
  write_ildg_gauge_conf(arg[6],conf,64,&mess);
  nissa_free(conf);
  
  ///////////////////////////////////////////

  closeNissa();

  return 0;
}
