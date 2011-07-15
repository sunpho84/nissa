#include <include.h>
#include <iostream>

using namespace std;

int icombo(int f,int icorr)
{return 8*f+icorr;}

int main()
{
  int T=16,njack=10;
  char basepath[1024]="/home/francesco/QCD/LAVORI/RHO/data/";
  
  FILE *fin=open_file("input","r");
  
  int nb;
  fscanf(fin,"%d",&nb);
  
  jvec M(nb,njack);
  for(int ib=0;ib<nb;ib++)
    {
      int b;
      fscanf(fin,"%d",&b);
      
      char path[1024];
      sprintf(path,"%s/%03d",basepath,b);
      int f=0;//,icorr=5;
      jvec a[8];
      for(int icorr=0;icorr<8;icorr++)
	{
	  a[icorr]=jvec_load(combine("%s/results.dat",path).c_str(),T,njack,icombo(f,icorr)).simmetrized(1);
	  a[icorr].print_to_file(combine("/tmp/%d_%03d",icorr,b).c_str());
	}
      cout<<b<<" "<<constant_fit(effective_mass(a[5]+a[6]),3,8)<<endl;
    }
  
  return 0;
}
