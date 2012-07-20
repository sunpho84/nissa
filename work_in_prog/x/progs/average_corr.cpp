#include <stdio.h>

#include "nissa.h"

#include "../src/types/types.h"
#include "../src/routines/read_and_write.h"

corr16 *unav_corr;

//initialize the program
void init_calc(int narg,char **arg)
{
  //Basic mpi initialization
  init_nissa();
  
  if(rank_tot>1) crash("only available in scalar");
  if(narg<4) crash("use %s file_in T L",arg[0]);
  
  int T=atoi(arg[2]);
  int L=atoi(arg[3]);
  
  //init the grid
  init_grid(T,L);
  
  //allocate correlator
  unav_corr=nissa_malloc("corr",loc_vol,corr16);
}

//close the program
void close_calc()
{
  nissa_free(unav_corr);
  
  close_nissa();
}

int main(int narg,char **arg)
{
  init_calc(narg,arg);
  
  read_corr16(unav_corr,arg[1]);
  
  FILE *f00=open_file("corr00_tau32-0_L24_T48_free.dat","w");
  FILE *f07=open_file("corr07_tau32-0_L24_T48_free.dat","w");
  FILE *f0305=open_file("corr0305_tau32-0_L24_T48_free.dat","w");
  FILE *f0406=open_file("corr0406_tau32-0_L24_T48_free.dat","w");
  
  coords x;
  for(x[0]=0;x[0]<glb_size[0]/2;x[0]++)
    for(x[1]=0;x[1]<glb_size[1]/2;x[1]++)
      for(x[2]=0;x[2]<glb_size[2]/2;x[2]++)
	for(x[3]=0;x[3]<glb_size[3]/2;x[3]++)
	  {
	    //compute distance
	    double d=0;
	    for(int mu=0;mu<4;mu++) d+=x[mu]*x[mu];
	    int i=glblx_of_coord(x);
	    
	    //average
	    corr16 av;
	    spinspin_put_to_zero(*((spinspin*)&av));
	    for(int iperm=0;iperm<6;iperm++)
	      for(int ipar=0;ipar<16;ipar++)
		{
		  int pel[6][4]={{0,1,2,3},{0,2,3,1},{0,3,1,2},{0,1,3,2},{0,3,2,1},{0,2,1,3}};
		  
		  coords c;
		  for(int mu=0;mu<4;mu++)
		    {
		      int p=(ipar&(1<<mu));
		      c[mu]=x[pel[iperm][mu]];
		      if(p) c[mu]=(glb_size[mu]-c[mu])%glb_size[mu];
		    }
		  spinspin_summassign(*((spinspin*)&av),*((spinspin*)(unav_corr+glblx_of_coord(c))));
		}
	    spinspin_prodassign_double(*((spinspin*)&av),1.0/(16*6));
	    
	    //if(d<=70)
	      {
		fprintf(f00,"%d %d %d %d %16.16le\n",x[1],x[2],x[3],x[0],av[5][0]);
		fprintf(f07,"%d %d %d %d %16.16le\n",x[1],x[2],x[3],x[0],av[0][0]);
		fprintf(f0305,"%d %d %d %d %16.16le\n",x[1],x[2],x[3],x[0],av[1][0]);
		fprintf(f0406,"%d %d %d %d %16.16le\n",x[1],x[2],x[3],x[0],av[6][0]);
	      }
	  }
  
  fclose(f00);
  fclose(f07);
  fclose(f0305);
  fclose(f0406);
  
  close_calc();
  
  return 0;
}
