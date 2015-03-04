#include <nissa.hpp>
#define N 21

using namespace nissa;

typedef complex ON_t[N];
comm_t lx_ON_comm,eo_ON_comm;

//initialize cpn simulation
void init_cpn(const char *path)
{
  open_input(path);
  
  //read geometry
  int L;
  read_str_int("L",&L);
  //init_grid(L,L);
  
  //set_lx_comm(lx_ON_comm,sizeof(ON_t));
}

THREADABLE_FUNCTION_0ARG(atry)
{
  GET_THREAD_ID();
  //printf("%d ciccicicicic\n",thread_id);
}
THREADABLE_FUNCTION_END

//close everything
void close_cpn()
{
  printf("Closing CPN\n");
}

void in_main(int narg,char **arg)
{
  //init
  if(narg<2) crash("Use %s input",arg[0]);
  // init_cpn(arg[1]);
  
  //atry();
  
  //close
  close_cpn();
}

int main(int narg,char **arg)
{
  init_nissa_threaded(narg,arg,in_main);
  
  close_nissa();
    
  return 0;
}
