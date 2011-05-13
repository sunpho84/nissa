#include "include.h"
#include "kl3_common.cpp"

jvec load_disconnected(char *base_path,int im,int ik,int r)
{
  char path_conn[1024];
  char path_disc[1024];
  sprintf(path_conn,"%s/oPVmuPo-sss_conf.1.dat",base_path);
  sprintf(path_disc,"%s/oPoVmuPo-sss_conf.1.dat",base_path);
  
  int a=1,b=1,c=0;
  return read_three_points(path_disc, a,im,a, b,ik, c,r,c,0,0)-read_three_points(path_conn,a,im,a,b,ik,b,r,c,0,0);
}

int main()
{
  njack=10;
  read_input();
  
  cout<<load_disconnected(base_path,0,0,0);

  return 0;
}
