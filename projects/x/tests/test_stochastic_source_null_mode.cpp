#include <math.h>

#include "nissa.hpp"
using namespace nissa;
using namespace std;

complex *source;
complex *spectre;

//initialize the program
void init_test(int narg,char **arg)
{
  //Basic mpi initialization
  init_nissa(narg,arg);
  
  //init the grid
  init_grid(8,8);
  
  //start loc rnd gen
  start_loc_rnd_gen(1);
  
  //allocate source and spectre
  source=nissa_malloc("source",loc_vol,complex);
  spectre=nissa_malloc("spectre",loc_vol,complex);
  
  //reset spectre
  memset(spectre,0,sizeof(complex)*loc_vol);
}

//close the program
void close_test()
{
  nissa_free(source);
  nissa_free(spectre);
  
  close_nissa();
}

//generate source with 0 method
void generate_source_0_method()
{
  double nm=0;
  
  //generate the source and compute null mode
  NISSA_LOC_VOL_LOOP(ivol)
    {
      source[ivol][0]=rnd_get_pm_one(&(loc_rnd_gen[ivol]));
      source[ivol][1]=0;
      nm+=source[ivol][0];
    }
  nm=glb_reduce_double(nm)/glb_vol;
  
  //subtract null mode
  NISSA_LOC_VOL_LOOP(ivol) source[ivol][0]-=nm;
}

//generate source with 1 method
void generate_source_1_method()
{
  //generate the source without null mode
  NISSA_LOC_VOL_LOOP(ivol)
    {
      if(ivol<glb_vol/2) source[ivol][0]=1;
      else source[ivol][0]=-1;
      source[ivol][1]=0;
    }
  
  //perform nsweep change
  int nsweep=10;
  for(int isweep=0;isweep<nsweep;isweep++)
    NISSA_LOC_VOL_LOOP(ivol)
      {
	int jvol=(int)rnd_get_unif(&glb_rnd_gen,0,glb_vol);
	if(ivol!=jvol)
	  {
	    int temp=source[ivol][0];
	    source[ivol][0]=source[jvol][0];
	    source[jvol][0]=temp;
	  }
      }
}

//summ the spectre to the already computed
void summ_the_spectre()
{
  //pass to momentum
  fft4d((complex*)source,(complex*)source,1,+1,0);
  
  //summ it
  NISSA_LOC_VOL_LOOP(ivol) spectre[ivol][0]+=squared_complex_norm(source[ivol]);
}

void write_file(const char *path,complex *data)
{
  //compute density and average
  FILE *fout=open_file(path,"w");
  fprintf(fout,"@s0 line type 0\n");
  fprintf(fout,"@    s0 symbol 2\n");
  fprintf(fout,"@    s0 symbol size 0.20000\n");
  fprintf(fout,"@    s0 symbol color 2\n");
  
  //loop over triangle
  coords x;
  for(x[0]=0;x[0]<=glb_size[0]/2;x[0]++)
    for(x[1]=0;x[1]<=x[0];x[1]++)
      for(x[2]=0;x[2]<=x[1];x[2]++)
	for(x[3]=0;x[3]<=x[2];x[3]++)
	  {
	    //compute distance
	    double dist=0;
	    for(int mu=0;mu<4;mu++) dist+=x[mu]*x[mu];
	    
	    double av=0;
	    for(int iperm=0;iperm<24;iperm++)
	      for(int ipar=0;ipar<16;ipar++)
		{
		  int pel[24][4]=
		    {
		      {0,1,2,3},{0,2,3,1},{0,3,1,2},{0,1,3,2},{0,3,2,1},{0,2,1,3},
		      {1,0,2,3},{1,2,3,0},{1,3,0,2},{1,0,3,2},{1,3,2,0},{1,2,0,3},
		      {2,1,0,3},{2,0,3,1},{2,3,1,0},{2,1,3,0},{2,3,0,1},{2,0,1,3},
		      {3,1,2,0},{3,2,0,1},{3,0,1,2},{3,1,0,2},{3,0,2,1},{3,2,1,0},
		    };
		  
		  coords c;
		  for(int mu=0;mu<4;mu++)
		    {
		      int p=(ipar&(1<<mu));
		      c[mu]=x[pel[iperm][mu]];
		      if(p) c[mu]=(glb_size[mu]-c[mu])%glb_size[mu];
		    }
		  av+=data[glblx_of_coord(c)][0];
		}
	    av/=24*16;
	    
	    fprintf(fout,"%lg %lg\n",dist,av);
	  }
  
  fclose(fout);
}


//study 0 method
void study_0_method()
{
  //average the spectral density between nsource
  int nsource=1000;
  for(int iso=0;iso<nsource;iso++)
    {
      printf("sourc %d/%d\n",iso+1,nsource);
      generate_source_1_method();
      summ_the_spectre();
    }
  NISSA_LOC_VOL_LOOP(ivol)
    spectre[ivol][0]/=nsource;
  
  //write spectral density
  write_file("spectral_density_1",spectre);
  
  //compute and write autocorr function
  fft4d((complex*)source,(complex*)spectre,1,-1,1);
  write_file("autocorr_function_1",source);
}

int main(int narg,char **arg)
{
  init_test(narg,arg);
  
  study_0_method();
  
  close_test();
  
  return 0;
}
