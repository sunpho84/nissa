#include "include.h"
#include "kl3_common.cpp"

jvec load_averaged_two_points(char *base_path,int im1,int im2,int ik1,int ik2,int r1,int r2)
{
  char path[1024];
  sprintf(path,"%s/oPPo-ss_conf.1.dat",base_path);
  
  int ri=0;
  
  //read the 2 two points with opposite momentum
  return (read_two_points(path,im1,im2,ik1,ik2,r1,r2,ri)+
	  read_two_points(path,im1,im2,ik2,ik1,r1,r2,ri))*0.5;
}

int main()
{
  njack=10;
  read_input();
  
  //load the neutral "up" flavour
  int r1=0,r2=r1;

  //load the pion
  int im1=0,im2=1;

  //loop over triangular impulses
  for(int ik1=0;ik1<nmoms;ik1++)
    for(int ik2=ik1;ik2<nmoms;ik2++)
      {
        //calculate q2
        double q2=3*pow((theta[ik1]-theta[ik2])*2*M_PI/L,2);
        //load the correlation functions, averaged      
        jvec c=load_averaged_two_points(base_path,im1,im2,ik1,ik2,r1,r2);
	jvec cs=c.simmetrized(1);
	jvec cmeff=effective_mass(cs);
	jack m=constant_fit(cmeff,13,23);
	cout<<q2<<" "<<m<<endl;
      }	

  return 0;
}
