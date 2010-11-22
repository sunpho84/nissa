#include <iostream>
#include <math.h>
int BIG=2<<29;
using namespace std;

int main()
{
  int L,T,N;

  cout<<"L? ";
  cin>>L;
  cout<<"T? ";
  cin>>T;
  cout<<"N? ";
  cin>>N;


  int V=L*L*L*T;
  int VL=V/N;

  int log2VL=0;
  do log2VL++;
  while ((2<<log2VL)<VL);

  int lista_fac[log2VL];

  int nfatt=0;
  int fatt=2;
  int TVL=VL;

  while(TVL>1)
    {
      int div=TVL/fatt;
      int res=TVL-div*fatt;
      if(res!=0) fatt++;
      else 
	{
	  TVL=div;
	  lista_fac[nfatt]=fatt;
	  nfatt++;
	}
    }

  int qtonfatt=1;
  for(int ifatt=1;ifatt<=nfatt;ifatt++) qtonfatt*=4;

  int minsurf=BIG;
  int P[4],mP[4];
  for(int ic=0;ic<qtonfatt;ic++)
    {
      P[0]=P[1]=P[2]=P[3]=1;

      int tic=ic;
      for(int ifatt=0;ifatt<nfatt;ifatt++)
	{
	  int i=tic%4;
	  tic/=4;
	  P[i]*=lista_fac[ifatt];
	}

      int intvol=1;
      for(int id=0;id<4;id++)
	if(P[id]>2) intvol*=P[id]-2;
	else intvol=0;
      
      int surf=VL-intvol;

      if(surf<minsurf && (L%P[0]==0) && (L%P[1]==0) && (L%P[2]==0) && (T%P[3]==0) )
	{
	  minsurf=surf;
	  for(int id=0;id<4;id++) mP[id]=P[id];
	}
    }

  cout<<"Vloc(x,y,z,t):\t"<<mP[0]<<" "<<mP[1]<<" "<<mP[2]<<" "<<mP[3]<<endl;
  cout<<"Nproc(x,y,z,t):\t"<<L/mP[0]<<" "<<L/mP[1]<<" "<<L/mP[2]<<" "<<T/mP[3]<<endl;

  return 0;
}
