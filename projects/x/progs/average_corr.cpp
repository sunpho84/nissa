#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "nissa.hpp"
using namespace std;

#include "../src/types/types.hpp"
#include "../src/routines/read_and_write.hpp"

corr16 *unav_corr;
corr16 *temp_corr;

//initialize the program
void init_calc(int narg,char **arg)
{
  //Basic mpi initialization
  init_nissa(narg,arg);
  
  if(nranks>1) crash("only available in scalar");
  if(narg<3) crash("use %s T L",arg[0]);
  
  int T=atoi(arg[1]);
  int L=atoi(arg[2]);
  
  //init the grid
  init_grid(T,L);
  
  //allocate correlator
  unav_corr=nissa_malloc("unav_corr",loc_vol,corr16);
  temp_corr=nissa_malloc("temp_corr",loc_vol,corr16);
}

//return the angle and the dist
void compute_dist2_angle(int &d2,double &angle,coords x)
{
  d2=angle=0;
  for(int mu=0;mu<4;mu++)
    {   
      d2+=x[mu]*x[mu];
      angle+=x[mu];
    }
  angle=acos(angle/sqrt(4*d2))*180/M_PI;
}

//close the program
void close_calc()
{
  nissa_free(unav_corr);
  nissa_free(temp_corr);
  
  close_nissa();
}

void print_ave(const char *outf,const char *suff)
{
  char path[1024];
  sprintf(path,"%s/corr00_tau32-0_L%2d_T%d_%s.dat",outf,glb_size[1],glb_size[0],suff);
  FILE *f00=open_file(path,"w");
  sprintf(path,"%s/corr07_tau32-0_L%2d_T%d_%s.dat",outf,glb_size[1],glb_size[0],suff);
  FILE *f07=open_file(path,"w");
  sprintf(path,"%s/corr0305_tau32-0_L%2d_T%d_%s.dat",outf,glb_size[1],glb_size[0],suff);
  FILE *f0305=open_file(path,"w");
  sprintf(path,"%s/corr0406_tau32-0_L%2d_T%d_%s.dat",outf,glb_size[1],glb_size[0],suff);
  FILE *f0406=open_file(path,"w");
  
  coords x;
  for(x[0]=0;x[0]<=glb_size[0]/2;x[0]++)
    for(x[1]=0;x[1]<=glb_size[1]/2;x[1]++)
      for(x[2]=0;x[2]<=glb_size[2]/2;x[2]++)
	for(x[3]=0;x[3]<=glb_size[3]/2;x[3]++)
	  {
	    int d2;
            double angle;
            compute_dist2_angle(d2,angle,x);
	    
	    //average
	    corr16 av;
	    int nperm,npar;
	    spinspin_put_to_zero(*((spinspin*)&av));
	    for(int ig=0;ig<16;ig++)
	      {
		if(ig==0||ig==5)
		  {
		    nperm=6;
		    npar=16;
		  }
		else
		  {
		    nperm=6;
		    npar=16;
		  }
		
		for(int iperm=0;iperm<nperm;iperm++)
		  for(int ipar=0;ipar<npar;ipar++)
		    {
		      int pel[6][4]={{0,1,2,3},{0,2,3,1},{0,3,1,2},{0,1,3,2},{0,3,2,1},{0,2,1,3}};
		      
		      coords c;
		      for(int mu=0;mu<4;mu++)
			{
			  int p=(ipar&(1<<mu));
			  c[mu]=x[pel[iperm][mu]];
			  if(p) c[mu]=(glb_size[mu]-c[mu])%glb_size[mu];
			}
		      complex_summassign(av[ig],unav_corr[glblx_of_coord(c)][ig]);
		    }
		complex_prodassign_double(av[ig],1.0/(nperm*npar));
	      }
	    
	    //filter 
	    if(d2<=70)
	      {
		double S0=-av[0][0];
		double P5=av[5][0];
		double V4=av[4][0],V1=av[1][0],V2=av[2][0],V3=av[3][0];
		double A4=av[5][0],A1=av[6][0],A2=av[7][0],A3=av[8][0];
		double Vi=V1+V2+V3;
		double Ai=A1+A2+A3;
		fprintf(f00,"%2d %2d %2d %2d %16.16le\n",x[1],x[2],x[3],x[0],P5);
		fprintf(f07,"%2d %2d %2d %2d %16.16le\n",x[1],x[2],x[3],x[0],S0);
		fprintf(f0305,"%2d %2d %2d %2d %16.16le %16.16le %16.16le %16.16le %16.16le %16.16le\n",x[1],x[2],x[3],x[0],V1,V2,V3,V4,Vi,Vi+V4);
		fprintf(f0406,"%2d %2d %2d %2d %16.16le %16.16le %16.16le %16.16le %16.16le %16.16le\n",x[1],x[2],x[3],x[0],A1,A2,A3,A4,Ai,Ai+A4);
	      }
	  }
  
  fclose(f00);
  fclose(f07);
  fclose(f0305);
  fclose(f0406);
}

void read_unave(corr16 *out,const char *name)
{
  char path[1024];
  sprintf(path,"../../raw_corrections/%2d/%s_corr",glb_size[1],name);
  
  read_corr16(out,path);
}

void write_unave(const char *name,corr16 *out)
{
  char path[1024];
  sprintf(path,"%s_corr",name);
  
  write_corr16(path,out,64);
}

void summ_with_coeff(double c)
{
  double_vector_summ_double_vector_prod_double((double*)unav_corr,(double*)unav_corr,(double*)temp_corr,c,sizeof(corr16)/sizeof(double)*loc_vol);
}

int main(int narg,char **arg)
{
  double mcrit=40.444;// is the AOKI calculation
  //double mcrit=41.6487;
  double mass_coef=mcrit/(16*M_PI*M_PI);
  init_calc(narg,arg);
    
  //prepare the tree coor
  vector_reset(unav_corr);
  read_unave(temp_corr,"tree");
  summ_with_coeff(3);
  print_ave("./","free");

  //prepare the first order corr
  vector_reset(unav_corr);
  read_unave(temp_corr,"mass");
  summ_with_coeff(2*4.0*mass_coef);
  read_unave(temp_corr,"self");
  summ_with_coeff(2*4.0);
  read_unave(temp_corr,"tad");
  summ_with_coeff(2*4.0);
  read_unave(temp_corr,"exch");
  summ_with_coeff(1*4.0);
  print_ave("./","first");
  
  close_calc();
  
  return 0;
}
