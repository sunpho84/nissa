#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "nissa.hpp"
using namespace std;

#include "../src/types/types.hpp"
#include "../src/routines/read_and_write.hpp"

const int plot_d2_max=70;
const int PP=0,SS=1,VV=2,AA=3;
char corr_name[4][8]={"00","07","0305","0406"};
const int corr_entr[4]={5,0,1,6};
const int corr_mult[4]={1,1,1,1};
const int recompute_first=0;

int *demo_of_loclx,*loclx_of_demo,*npoints_dist2;
int *dist2;
int ndemo_points,dist2_max;
double beta,mass,g2;
const char base_path[]="/Users/francesco/QCD/LAVORI/X";

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
void prepare_demo_table(double cut_angle,int max_dist)
{
  NISSA_LOC_VOL_LOOP(ivol) demo_of_loclx[ivol]=-1;
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
	    if(angle<=(cut_angle+1.e-6) && check_not_too_far(x,max_dist))
	      {
		demo_of_loclx[ivol]=ndemo_points++;
		
		//compute point degeneracy
		npoints_dist2[dist2[ivol]]++;
	      }
	  }
  
  //allocate back-mapping table and fill
  loclx_of_demo=nissa_malloc("loclx_of_demo",ndemo_points,int);
  NISSA_LOC_VOL_LOOP(ivol)
    if(demo_of_loclx[ivol]!=-1) loclx_of_demo[demo_of_loclx[ivol]]=ivol;
  
  master_printf("Number of democratic points: %d\n",ndemo_points);
}

//initialize the program
void init_calc(int narg,char **arg)
{
  //Basic mpi initialization
  init_nissa(narg,arg);
  
  if(nranks>1) crash("only available in scalar");
  if(narg<2) crash("use %s input",arg[0]);
  
  open_input(arg[1]);
  int L;
  read_str_int("L",&L);
  int T;
  read_str_int("T",&T);
    
  double cut_angle;
  read_str_double("CutAngle",&cut_angle);
  
  int cut_dist;
  read_str_int("CutDist",&cut_dist);
  
  read_str_double("Beta",&beta);
  read_str_double("Mass",&mass);
  
  double plaq;
  read_str_double("Plaq",&plaq);
  g2=6/beta/plaq;
  master_printf("g2_boost: %lg\n",g2);

  //init the grid
  init_grid(T,L);
  
  //allocate table of distances
  dist2=nissa_malloc("dist2",loc_vol,int);
  
  //allocate the table of degeneracy of points
  dist2_max=compute_dist2_max(); 
  npoints_dist2=nissa_malloc("npoints_dist2",dist2_max,int);
  
  //allocate table of democratic points
  demo_of_loclx=nissa_malloc("demo_of_loclcx",loc_vol,int);
  prepare_demo_table(cut_angle,cut_dist);
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

double Zp_minus_one(double d2,int icorr)
{
  double Cf=4.0/3;
  double tilde2=Cf/(16*M_PI*M_PI);
  double egamma=0.5772156649015329;
  
  double Zpmo;
  
  switch(icorr)
    {
    case PP:
      Zpmo=(-18.7333556530313-3*log(d2/4.0)-6*egamma);
      break;
    case SS:
      Zpmo=(-10.606730862152009-3*log(d2/4.0)-6*egamma);
      break;
    case VV:
      Zpmo=-18.102886492137763;
      break;
    case AA:
      Zpmo=-14.039574096698127;
      break;
    }
  
  return tilde2*Zpmo;
}


//load the real part of democratic points of the correlator
void load_demo_averaged_text_corr(double *out,char *path)
{
  master_printf("Opening: %d\n",path);
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

  int nperm=6;
  int npar=16;
  
  if(average) master_printf("Averaging\n");
  
  //copy only democratic points
  NISSA_LOC_VOL_LOOP(ivol)
  {
    int idemo=demo_of_loclx[ivol];
    if(idemo!=-1)
      {
	if(average)
	  {
	    spinspin av;
	    spinspin_put_to_zero(av);
	    //takes all the copies
	    for(int iperm=0;iperm<nperm;iperm++)
	      for(int ipar=0;ipar<npar;ipar++)
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
	    spinspin_prod_double(*((spinspin*)&out[idemo]),*((spinspin*)&av),1.0/(nperm*npar));
	  }
	else memcpy(out[idemo],temp[ivol],sizeof(corr16));
	
	//put 1+2+3+4 in 1
	for(int imu=2;imu<=4;imu++)
	  complex_summassign(out[idemo][1],out[idemo][imu]);
	//put 6+7+8+9 in 6
	for(int imu=7;imu<=9;imu++)
	  complex_summassign(out[idemo][6],out[idemo][imu]);
      }
  }
  
  nissa_free(temp);
}

void load_correction(corr16 *out,const char *inf,int icorr,const char *suff)
{
  vector_reset(out);
  
  char path[1024];
  sprintf(path,"%s/corr%s_tau32-0_L%2d_T%2d_%s.dat",inf,corr_name[icorr],glb_size[1],glb_size[0],suff);
  master_printf("Opening: %s\n",path);
  FILE *fin=open_file(path,"r");
  
  coords x;
  for(x[0]=0;x[0]<=glb_size[0]/2;x[0]++)
    for(x[1]=0;x[1]<=glb_size[1]/2;x[1]++)
      for(x[2]=0;x[2]<=glb_size[2]/2;x[2]++)
        for(x[3]=0;x[3]<=glb_size[3]/2;x[3]++)
	  {
	    int ivol=glblx_of_coord(x);
	    int idemo=demo_of_loclx[ivol];

	    int d2;
	    double angle;
	    compute_dist2_angle(d2,angle,x);
	    
	    //read only what really wrote
	    if(d2<=70) //only 70 present
	      {
		coords y;
		double temp;
		switch(icorr)
		  {
		    double dum1,dum2,dum3,dum4,dum5;
		    int nr;
		  case PP:
		    nr=fscanf(fin,"%d %d %d %d %lg",&y[1],&y[2],&y[3],&y[0],&temp);
		    if(nr!=5) crash("read %d",nr);
		    break;
		  case SS:
		    nr=fscanf(fin,"%d %d %d %d %lg",&y[1],&y[2],&y[3],&y[0],&temp);
		    if(nr!=5) crash("read %d",nr);
		    break;
		  case VV:
		    nr=fscanf(fin,"%d %d %d %d %lg %lg %lg %lg %lg %lg",&y[1],&y[2],&y[3],&y[0],&dum1,&dum2,&dum3,&dum4,&dum5,&temp);
		    if(nr!=10) crash("read %d",nr);
		    break;
		  case AA:
		    nr=fscanf(fin,"%d %d %d %d %lg %lg %lg %lg %lg %lg",&y[1],&y[2],&y[3],&y[0],&dum1,&dum2,&dum3,&dum4,&dum5,&temp);
		    if(nr!=10) crash("read %d",nr);
		    break;
		  }
		
		if(idemo!=-1)
		  out[idemo][corr_entr[icorr]][0]=temp;
	      }
	  }
  
  char fuf[1024];
  if(fscanf(fin,"%s",fuf)==1) crash("not reached EOF!");
  
  fclose(fin);
}

//write parameters of plot
void write_pars(FILE *fout)
{
  fprintf(fout,"@s0 line type 0\n@s0 symbol 1\n@s0 symbol size 0.240000\n@s0 symbol fill pattern 1\n");
  fprintf(fout,"@world 0, 1e-07, %d, 0.1\n@yaxes scale Logarithmic\n@yaxis  tick major 10\n@yaxis  tick minor ticks 9\n",plot_d2_max);
}

//write directly not averaging
void write_demo(const char *path,double *c)
{
  FILE *fout=open_file(path,"w");
  
  //write parameters
  write_pars(fout);
  
  for(int idemo=0;idemo<ndemo_points;idemo++)
    {
      int ivol=loclx_of_demo[idemo];
      if(dist2[ivol]<=plot_d2_max)
	if(glb_coord_of_loclx[ivol][3]>=glb_coord_of_loclx[ivol][2]&&glb_coord_of_loclx[ivol][2]>=glb_coord_of_loclx[ivol][1])
	  if(!isnan(c[idemo]))
	    fprintf(fout,"%lg %lg\n",(double)(dist2[ivol]),c[idemo]);
    }

  fclose(fout);
}

//write averaging points with the same distance
void write_allaveraged_demo(const char *path,double *c)
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
  for(int d2=0;d2<=plot_d2_max;d2++)
    if(npoints_dist2[d2]!=0)
      if(!isnan(cave[d2]))
	fprintf(fout,"%d %lg\n",d2,cave[d2]);  

  fclose(fout);
}

void correct(int icorr)
{
  //load uncorrected data
  double *full=nissa_malloc("full",ndemo_points,double);
  char uncorr_path[1024];
  sprintf(uncorr_path,"%s/uncorrected_data/%1.2f/%2d/%1.4f/corr%s_tau32-0_b%1.2f_mu%1.4f_L%2d_T%2d.dat",base_path,beta,glb_size[1],mass,corr_name[icorr],beta,mass,glb_size[1],glb_size[0]);
  load_demo_averaged_text_corr(full,uncorr_path);
  
  //load the tree level correlation on lattice
  char path[1024];
  corr16 *tree_corr_lat=nissa_malloc("tree_corr_lat",ndemo_points,corr16);
  sprintf(path,"%s/corrections/%d/",base_path,glb_size[1]);
  load_correction(tree_corr_lat,path,icorr,"free");

  sprintf(path,"%s/raw_corrections/%d/tree_corr",base_path,glb_size[1]);
  if(recompute_first) load_demo_ildg_corr(tree_corr_lat,path);
  
  //load the first order correlation on lattice
  corr16 *first_corr_lat=nissa_malloc("first_corr_lat",ndemo_points,corr16);
  sprintf(path,"%s/corrections/%d/",base_path,glb_size[1]);
  load_correction(first_corr_lat,path,icorr,"first");
  
  corr16 *mass_corr_lat=nissa_malloc("mass_corr_lat",ndemo_points,corr16);
  sprintf(path,"%s/raw_corrections/%d/mass_corr",base_path,glb_size[1]);
  if(recompute_first) load_demo_ildg_corr(mass_corr_lat,path,true);
  
  corr16 *exch_corr_lat=nissa_malloc("exch_corr_lat",ndemo_points,corr16);
  sprintf(path,"%s/raw_corrections/%d/exch_corr",base_path,glb_size[1]);
  if(recompute_first) load_demo_ildg_corr(exch_corr_lat,path,true);
  
  corr16 *self_corr_lat=nissa_malloc("self_corr_lat",ndemo_points,corr16);
  sprintf(path,"%s/raw_corrections/%d/self_corr",base_path,glb_size[1]);
  if(recompute_first) load_demo_ildg_corr(self_corr_lat,path);

  corr16 *tad_corr_lat=nissa_malloc("tad_corr_lat",ndemo_points,corr16);
  sprintf(path,"%s/raw_corrections/%d/tad_corr",base_path,glb_size[1]);
  if(recompute_first) load_demo_ildg_corr(tad_corr_lat,path);
  
  //correlations
  double *tree_cont=nissa_malloc("tree_cont",ndemo_points,double);
  double *tree_lat=nissa_malloc("tree_lat",ndemo_points,double);
  double *slope_cont=nissa_malloc("slope_cont",ndemo_points,double);
  double *slope_lat=nissa_malloc("slope_lat",ndemo_points,double);
  double *tree_diff=nissa_malloc("tree_diff",ndemo_points,double);
  double *tree_ratio=nissa_malloc("tree_ratio",ndemo_points,double);
  double *full_tree_subt=nissa_malloc("tree_subt",ndemo_points,double);
  double *full_tree_divided=nissa_malloc("tree_divided",ndemo_points,double);
  double *first_in_lat=nissa_malloc("first_lat",ndemo_points,double);
  double *first_lat=nissa_malloc("first_lat",ndemo_points,double);
  double *mass_lat=nissa_malloc("mass_lat",ndemo_points,double);
  double *exch_lat=nissa_malloc("exch_lat",ndemo_points,double);
  double *tad_lat=nissa_malloc("tad_lat",ndemo_points,double);
  double *self_lat=nissa_malloc("self_lat",ndemo_points,double);
  double *first_cont=nissa_malloc("first_cont",ndemo_points,double);
  double *first_diff=nissa_malloc("first_diff",ndemo_points,double);
  double *first_ratio=nissa_malloc("first_ratio",ndemo_points,double);
  double *tree_plus_first_lat=nissa_malloc("tree_plus_first_lat",ndemo_points,double);
  double *tree_plus_first_cont=nissa_malloc("tree_plus_first_cont",ndemo_points,double);
  double *full_first_subt=nissa_malloc("full_first_subt",ndemo_points,double);
  double *full_first_divided=nissa_malloc("full_first_divided",ndemo_points,double);
  
  double *Z_uncorr=nissa_malloc("Z_uncorr",ndemo_points,double);
  double *Z_tree_subt=nissa_malloc("Z_tree_subt",ndemo_points,double);
  double *Z_tree_divided=nissa_malloc("Z_tree_divided",ndemo_points,double);
  double *Z_first_divided=nissa_malloc("Z_first_divided",ndemo_points,double);
  
  for(int idemo=0;idemo<ndemo_points;idemo++)
    {
      //add sign to full
      if(icorr==SS) full[idemo]*=-1;
      
      //take dist2
      int ivol=loclx_of_demo[idemo];
      int d2=dist2[ivol];
      
      //compute tree continuum function
      const int cont_coeff[4]={3,3,6,6};
      const int sign_corr[4]={1,1,1,1};
      tree_cont[idemo]=cont_coeff[icorr]/(M_PI*M_PI*M_PI*M_PI*d2*d2*d2);
      
      //compute Z uncorr
      Z_uncorr[idemo]=sqrt(tree_cont[idemo]/full[idemo]);
      
      ////////////////////////////////////// tree corr //////////////////////////////////////////
      
      //compute tree lattice
      tree_lat[idemo]=sign_corr[icorr]*tree_corr_lat[idemo][corr_entr[icorr]][RE];
      if(recompute_first) tree_lat[idemo]=sign_corr[icorr]*cont_coeff[icorr]*tree_corr_lat[idemo][corr_entr[icorr]][RE];

      //diff between tree lat and tree cont
      tree_diff[idemo]=tree_lat[idemo]-tree_cont[idemo];
      
      //ratio between tree lat and tree cont
      tree_ratio[idemo]=tree_lat[idemo]/tree_cont[idemo];

      //subtract tree level discretization
      full_tree_subt[idemo]=full[idemo]-pow(Z_uncorr[icorr],2)*tree_diff[idemo];
      
      //divide by tree level discretization-continuum ratio
      full_tree_divided[idemo]=full[idemo]/tree_ratio[idemo];
      
      //compute subtracted Z
      Z_tree_subt[idemo]=sqrt(tree_cont[idemo]/full_tree_subt[idemo]);
      
      //compute divided Z
      Z_tree_divided[idemo]=sqrt(tree_cont[idemo]/full_tree_divided[idemo]);
            
      ////////////////////////////////////////// first //////////////////////////////////////////
      
      //compute first order in the continuum
      first_cont[idemo]=-2*tree_cont[idemo]*Zp_minus_one(d2,icorr);
      
      //load corrections
      mass_lat[idemo]=mass_corr_lat[idemo][corr_entr[icorr]][RE];
      self_lat[idemo]=self_corr_lat[idemo][corr_entr[icorr]][RE];
      tad_lat[idemo]=tad_corr_lat[idemo][corr_entr[icorr]][RE];
      exch_lat[idemo]=exch_corr_lat[idemo][corr_entr[icorr]][RE];
      
      //slope
      slope_lat[idemo]=self_lat[idemo]/tree_lat[idemo];
      slope_cont[idemo]=1+g2*Zp_minus_one(d2,icorr);
      
      //compute first order on lattice
      double mass_coef=0.25611505461400043720;
      first_lat[idemo]=sign_corr[icorr]*first_corr_lat[idemo][corr_entr[icorr]][RE];
      if(recompute_first) first_lat[idemo]=sign_corr[icorr]*4*(exch_lat[idemo]+2*self_lat[idemo]+2*tad_lat[idemo]+2*mass_coef*mass_lat[idemo]);
      //if(recompute_first) first_lat[idemo]=sign_corr[icorr]*4*(exch_lat[idemo]+2*self_lat[idemo]-2*tad_lat[idemo]);
      
      //diff between first lat and first cont
      first_diff[idemo]=first_lat[idemo]-first_cont[idemo];      
      
      //diff between first lat and first cont
      full_first_subt[idemo]=full[idemo]-first_diff[idemo];
      
      //tree_plus_first_lat
      tree_plus_first_lat[idemo]=tree_lat[idemo]+g2*first_lat[idemo];
      tree_plus_first_cont[idemo]=tree_cont[idemo]+g2*first_cont[idemo];
      
      //ratio between tree+first lat and tree+first cont
      first_ratio[idemo]=tree_plus_first_lat[idemo]/tree_plus_first_cont[idemo];
      
      //correct with ratio
      full_first_divided[idemo]=full[idemo]/first_ratio[idemo];
      
      //compute divided Z
      Z_first_divided[idemo]=sqrt(tree_cont[idemo]/full_first_divided[idemo]);
    }
  
  write_demo("plots/full.xmg",full);
  write_allaveraged_demo("plots/tree_cont.xmg",tree_cont);
  write_demo("plots/tree_lat.xmg",tree_lat);
  write_demo("plots/tree_diff.xmg",tree_diff);
  write_demo("plots/tree_ratio.xmg",tree_ratio);
  write_demo("plots/full_tree_subtracted.xmg",full_tree_subt);
  write_demo("plots/full_tree_divided.xmg",full_tree_divided);

  write_allaveraged_demo("plots/first_cont.xmg",first_cont);
  write_demo("plots/slope_lat.xmg",slope_lat);
  write_demo("plots/slope_cont.xmg",slope_cont);
  write_demo("plots/mass_lat.xmg",mass_lat);
  write_demo("plots/self_lat.xmg",self_lat);
  write_demo("plots/tad_lat.xmg",tad_lat);
  write_demo("plots/exch_lat.xmg",exch_lat);
  write_demo("plots/first_in_lat.xmg",first_in_lat);
  write_demo("plots/first_lat.xmg",first_lat);
  write_demo("plots/first_diff.xmg",first_diff);
  write_demo("plots/first_ratio.xmg",first_ratio);
  write_demo("plots/tree_plus_first_lat.xmg",tree_plus_first_lat);
  write_demo("plots/tree_plus_first_cont.xmg",tree_plus_first_cont);
  write_demo("plots/full_first_subtracted.xmg",full_first_subt);
  write_demo("plots/full_first_divided.xmg",full_first_divided);
  
  write_demo("plots/Z_uncorr.xmg",Z_uncorr);
  write_demo("plots/Z_tree_subt.xmg",Z_tree_subt);
  write_demo("plots/Z_tree_divided.xmg",Z_tree_divided);
  write_demo("plots/Z_first_divided.xmg",Z_first_divided);

  write_allaveraged_demo("plots/Z_uncorr_ave.xmg",Z_uncorr);
  write_allaveraged_demo("plots/Z_tree_subt_ave.xmg",Z_tree_subt);
  write_allaveraged_demo("plots/Z_tree_divided_ave.xmg",Z_tree_divided);
  write_allaveraged_demo("plots/Z_first_divided_ave.xmg",Z_first_divided);
  
  nissa_free(full);
  nissa_free(tree_corr_lat);
  nissa_free(tree_cont);
  nissa_free(tree_lat);
  nissa_free(mass_lat);
  nissa_free(slope_lat);
  nissa_free(slope_cont);
  nissa_free(tree_diff);
  nissa_free(tree_ratio);
  nissa_free(full_tree_subt);
  nissa_free(full_tree_divided);
  nissa_free(first_corr_lat);
  nissa_free(mass_corr_lat);
  nissa_free(exch_corr_lat);
  nissa_free(tad_corr_lat);
  nissa_free(self_corr_lat);
  nissa_free(first_lat);
  nissa_free(first_in_lat);
  nissa_free(self_lat);
  nissa_free(exch_lat);
  nissa_free(tad_lat);
  nissa_free(first_cont);
  nissa_free(first_diff);
  nissa_free(first_ratio);
  nissa_free(tree_plus_first_lat);
  nissa_free(tree_plus_first_cont);
  nissa_free(full_first_subt);
  nissa_free(full_first_divided);
  
  nissa_free(Z_uncorr);
  nissa_free(Z_tree_subt);
  nissa_free(Z_tree_divided);
  nissa_free(Z_first_divided);
}

int main(int narg,char **arg)
{
  init_calc(narg,arg);
  
  correct(PP);
  //correct(SS);
  //correct(VV);
  //correct(AA);
  
  close_calc();
  
  return 0;
}
