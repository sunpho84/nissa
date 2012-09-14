#include <stdio.h>
#include <math.h>

#include "nissa.h"

#include "../src/types/types.h"
#include "../src/routines/read_and_write.h"

char P5P5_filename[]="corr00";
int *demo_of_loclx,*loclx_of_demo,*npoints_dist2;
int *dist2;
int ndemo_points,dist2_max;

//compute the maximal distance
int compute_dist2_max()
{
  int dist2_max=0;
  for(int mu=0;mu<4;mu++) dist2_max+=sqr(glb_size[mu]/2);
  
  return dist2_max+1;
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

//return whether the point is not too far
int check_not_too_far(coords x,int max_dist)
{
  int not_too_far=1;
  for(int mu=0;mu<4;mu++)
    if(x[mu]>max_dist) not_too_far=0;
  
  return not_too_far;
}

//select democratic points
void prepare_demo_table(double cut_angle,double max_dist)
{
  nissa_loc_vol_loop(ivol) demo_of_loclx[ivol]=-1;
  vector_reset(npoints_dist2);
  
  //loop on triangle
  coords x;
  ndemo_points=0;
  for(x[0]=0;x[0]<=glb_size[0]/2;x[0]++)
    for(x[1]=0;x[1]<=glb_size[1]/2;x[1]++)
      for(x[2]=x[1];x[2]<=glb_size[2]/2;x[2]++)
	for(x[3]=x[2];x[3]<=glb_size[3]/2;x[3]++)
	  {
	    //take index of coords
	    int ivol=glblx_of_coord(x);
	    
	    //compute distance and angle, and check if not too far
	    double angle;
	    compute_dist2_angle(dist2[ivol],angle,x);
	    
	    //check if democratic and not too far
	    if(angle<=cut_angle && check_not_too_far(x,max_dist))
	      {
		demo_of_loclx[ivol]=ndemo_points++;
		
		//compute point degeneracy
		npoints_dist2[dist2[ivol]]++;
	      }
	  }
  
  //allocate back-mapping table and fill
  loclx_of_demo=nissa_malloc("loclx_of_demo",ndemo_points,int);
  nissa_loc_vol_loop(ivol)
    if(demo_of_loclx[ivol]!=-1) loclx_of_demo[demo_of_loclx[ivol]]=ivol;
  
  master_printf("Number of democratic points: %d\n",ndemo_points);
}

//initialize the program
void init_calc(int narg,char **arg)
{
  //Basic mpi initialization
  init_nissa();
  
  if(nissa_nranks>1) crash("only available in scalar");
  if(narg<5) crash("use %s T L base angle",arg[0]);
  
  int T=atoi(arg[1]);
  int L=atoi(arg[2]);
  
  double cut_angle;
  sscanf(arg[4],"%lg",&cut_angle);
  
  //init the grid
  init_grid(T,L);
  
  //allocate table of distances
  dist2=nissa_malloc("dist2",loc_vol,int);
  
  //allocate the table of degeneracy of points
  dist2_max=compute_dist2_max(); 
  npoints_dist2=nissa_malloc("npoints_dist2",dist2_max,int);
  
  //allocate table of democratic points
  demo_of_loclx=nissa_malloc("demo_of_loclcx",loc_vol,int);
  prepare_demo_table(cut_angle,40);
}

//close the program
void close_calc()
{
  nissa_free(dist2);
  nissa_free(npoints_dist2);
  nissa_free(demo_of_loclx);
  nissa_free(loclx_of_demo);
  
  close_nissa();
}

//load the real part of democratic points of the correlator
void load_demo_averaged_text_corr(double *out,char *path)
{
  FILE *fin=open_file(path,"r");
  
  coords x;
  for(x[0]=0;x[0]<=glb_size[0]/2;x[0]++)
    for(x[1]=0;x[1]<=glb_size[1]/2;x[1]++)
      for(x[2]=x[1];x[2]<=glb_size[2]/2;x[2]++)
	for(x[3]=x[2];x[3]<=glb_size[3]/2;x[3]++)
	  {
	    int ivol=glblx_of_coord(x);
	    double angle;
	    int d2;
	    compute_dist2_angle(d2,angle,x);
	    
	    //load
	    coords y;
	    double yd2,yangle;
	    double re,im;
	    int nr=fscanf(fin,"%d %d %d %d %lg %lg %lg %lg",&y[1],&y[2],&y[3],&y[0],&yd2,&yangle,&re,&im);
	    
	    //checks
	    if(nr!=8) crash("Read %d instead than 8",nr);
	    for(int mu=0;mu<4;mu++) if(y[mu]!=x[mu]) crash("Read %d instead than %d for dir %d",y[mu],x[mu],mu);
	    if(int(sqrt(yd2))!=int(sqrt(dist2[ivol]))) crash("Distance read %lg does not coincide with expected %lg",yd2,dist2[ivol]);
	    if(!isnan(angle)&&fabs(angle-yangle)>1.e-4) crash("Angle expected: %lg, read: %lg",angle,yangle);
	    
	    //check democracy
	    int idemo=demo_of_loclx[ivol];
	    if(idemo!=-1)
	      out[idemo]=re;
	  }
  
  fclose(fin);
}

//load the democratic points of an ildg correlator
void load_demo_ildg_corr(corr16 *out,char *path,bool average=false)
{
  //load all points
  corr16 *temp=nissa_malloc("temp",loc_vol,corr16);
  read_corr16(temp,path);

  //copy only democratic points
  nissa_loc_vol_loop(ivol)
    {
      int idemo=demo_of_loclx[ivol];
      if(idemo!=-1)
	if(average)
	  {
	    spinspin av;
	    spinspin_put_to_zero(av);
	    //takes all the copies
	    for(int iperm=0;iperm<6;iperm++)
	      for(int ipar=0;ipar<16;ipar++)
		{
		  int pel[6][4]={{0,1,2,3},{0,2,3,1},{0,3,1,2},{0,1,3,2},{0,3,2,1},{0,2,1,3}};
		  
		  coords c;
		  for(int mu=0;mu<4;mu++)
		    {
		      int p=(ipar&(1<<mu));
		      c[mu]=glb_coord_of_loclx[ivol][pel[iperm][mu]];
		      if(p) c[mu]=(glb_size[mu]-c[mu])%glb_size[mu];
		    }
		  spinspin_summassign(*((spinspin*)&av),*((spinspin*)(temp+glblx_of_coord(c))));
		}
	    spinspin_prod_double(*((spinspin*)&out[idemo]),*((spinspin*)&av),1.0/(16*6));
	  }
      else memcpy(out[idemo],temp[ivol],sizeof(corr16));
    }
  
  if(average)
    {
    }
  
  nissa_free(temp);
}

//write parameters of plot
void write_pars(FILE *fout)
{
  fprintf(fout,"@s0 line type 0\n@s0 symbol 1\n@s0 symbol size 0.240000\n@s0 symbol fill pattern 1\n");
  fprintf(fout,"@world 0, 1e-07, 70, 0.1\n@yaxes scale Logarithmic\n@yaxis  tick major 10\n@yaxis  tick minor ticks 9\n");
}

//write directly not averaging
void write_unaveraged_demo(const char *path,double *c)
{
  FILE *fout=open_file(path,"w");
  
  //write parameters
  write_pars(fout);
  
  for(int idemo=0;idemo<ndemo_points;idemo++)
    fprintf(fout,"%d %lg\n",dist2[loclx_of_demo[idemo]],c[idemo]);

  fclose(fout);
}

//write averaging points with the same distance
void write_averaged_demo(const char *path,double *c)
{
  FILE *fout=open_file(path,"w");
  
  //write parameters
  write_pars(fout);
  
  //reset the average
  double cave[dist2_max];
  memset(cave,0,sizeof(double)*dist2_max);
  
  //average
  for(int idemo=0;idemo<ndemo_points;idemo++)
    {
      int ivol=loclx_of_demo[idemo];
      int d2=dist2[ivol];
      cave[d2]+=c[idemo]/npoints_dist2[d2];
    }
  
  //print
  for(int d2=0;d2<dist2_max;d2++)
    if(npoints_dist2[d2]!=0)
      fprintf(fout,"%d %lg\n",d2,cave[d2]);  

  fclose(fout);
}

int main(int narg,char **arg)
{
  init_calc(narg,arg);

  //load uncorrected data
  double *full=nissa_malloc("full",ndemo_points,double);;
  char uncorr_path[1024];
  sprintf(uncorr_path,"correlations/%s%s",P5P5_filename,arg[3]);
  load_demo_averaged_text_corr(full,uncorr_path);
  
  //load the tree level diagram
  corr16 *tree_diag=nissa_malloc("tree_diagr",ndemo_points,corr16);
  load_demo_ildg_corr(tree_diag,(char*)"corrections/tree_corr");
  
  //load the self energy diagram
  corr16 *self_diag=nissa_malloc("self_diag",ndemo_points,corr16);
  load_demo_ildg_corr(self_diag,(char*)"corrections/self_corr");
  
  //load the tadpole diagram
  corr16 *tad_diag=nissa_malloc("tad_diag",ndemo_points,corr16);
  load_demo_ildg_corr(tad_diag,(char*)"corrections/tad_corr");
  
  //load the exchange diagram
  corr16 *exch_diag=nissa_malloc("exch_diag",ndemo_points,corr16);
  load_demo_ildg_corr(exch_diag,(char*)"corrections/exch_corr",true);
  
  //correlations
  double *tree_cont=nissa_malloc("tree_cont",ndemo_points,double);
  double *tree_lat=nissa_malloc("tree_lat",ndemo_points,double);
  double *tree_diff=nissa_malloc("tree_diff",ndemo_points,double);
  double *tree_ratio=nissa_malloc("tree_ratio",ndemo_points,double);
  double *full_tree_subt=nissa_malloc("tree_subt",ndemo_points,double);
  double *full_tree_divided=nissa_malloc("tree_divided",ndemo_points,double);
  double *self_lat=nissa_malloc("self_lat",ndemo_points,double);
  double *tad_lat=nissa_malloc("tad_lat",ndemo_points,double);
  double *exch_lat=nissa_malloc("exch_lat",ndemo_points,double);
  double *first_lat=nissa_malloc("first_lat",ndemo_points,double);
  double *first_cont=nissa_malloc("first_cont",ndemo_points,double);
  double *first_diff=nissa_malloc("first_diff",ndemo_points,double);
  
  for(int idemo=0;idemo<ndemo_points;idemo++)
    {
      //take dist2
      int ivol=loclx_of_demo[idemo];
      int d2=dist2[ivol];
      
      //compute tree continuum function
      tree_cont[idemo]=3/(M_PI*M_PI*M_PI*M_PI*d2*d2*d2);
      
      //compute tree lattice
      tree_lat[idemo]=tree_diag[idemo][5][RE]*3;
      
      //diff between tree lat and tree cont
      tree_diff[idemo]=tree_lat[idemo]-tree_cont[idemo];
      
      //ratio between tree lat and tree cont
      tree_ratio[idemo]=tree_lat[idemo]/tree_cont[idemo];

      //subtract tree level discretization
      full_tree_subt[idemo]=full[idemo]-tree_diff[idemo];
      
      //divide by tree level discretization-continuum ratio
      full_tree_divided[idemo]=full[idemo]/tree_ratio[idemo];
      
      ////////////// first ///////////
      
      //constants
      double tilde2=1/(16*M_PI*M_PI);
      double egamma=0.5772156649015329;
      ///////////////////// attention: changed sign to log
      double Zp_minus_one=tilde2*(-18.7334+3*log(d2/4.0)-6*egamma);
      
      //compute first order in the continuum
      first_cont[idemo]=-tree_cont[idemo]*2*Zp_minus_one;
      
      //compute self, tad and exch
      self_lat[idemo]=self_diag[idemo][5][RE]*3;
      tad_lat[idemo]=tad_diag[idemo][5][RE]*3;
      exch_lat[idemo]=exch_diag[idemo][5][RE]*3;
      
      //compute first order on lattice
      first_lat[idemo]=0;
      first_lat[idemo]-=self_diag[idemo][5][RE]*4*2;
      first_lat[idemo]-=tad_diag[idemo][5][RE]*4*2;
      first_lat[idemo]-=exch_diag[idemo][5][RE]*4;
      
      //diff between first lat and first cont
      first_diff[idemo]=first_lat[idemo]-first_cont[idemo];      
    }
  
  write_unaveraged_demo("plots/full.xmg",full);
  write_unaveraged_demo("plots/tree_cont.xmg",tree_cont);
  write_unaveraged_demo("plots/tree_lat.xmg",tree_lat);
  write_unaveraged_demo("plots/tree_diff.xmg",tree_diff);
  write_unaveraged_demo("plots/tree_ratio.xmg",tree_ratio);
  write_unaveraged_demo("plots/full_tree_subtracted.xmg",full_tree_subt);
  write_unaveraged_demo("plots/full_tree_divided.xmg",full_tree_divided);
  
  write_unaveraged_demo("plots/first_cont.xmg",first_cont);
  write_unaveraged_demo("plots/self_lat.xmg",self_lat);
  write_unaveraged_demo("plots/tad_lat.xmg",tad_lat);
  write_unaveraged_demo("plots/exch_lat.xmg",exch_lat);
  write_unaveraged_demo("plots/first_lat.xmg",first_lat);
  write_unaveraged_demo("plots/first_diff.xmg",first_diff);
  
  nissa_free(full);
  nissa_free(tree_diag);
  nissa_free(self_diag);
  nissa_free(tad_diag);
  nissa_free(exch_diag);
  nissa_free(self_lat);
  nissa_free(tad_lat);
  nissa_free(exch_lat);
  nissa_free(tree_cont);
  nissa_free(tree_lat);
  nissa_free(tree_diff);
  nissa_free(tree_ratio);
  nissa_free(full_tree_subt);
  nissa_free(full_tree_divided);
  nissa_free(first_lat);
  nissa_free(first_cont);
  nissa_free(first_diff);
  
  close_calc();
  //1+g2tilde*(-18.7334-3log(x2/4)-6gammaE)
  return 0;
}
