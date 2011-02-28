#pragma once

//search the index of the theta with opposite sign
int opposite_theta(int ik,double *theta)
{
  int ik_opp=0;
  
  while(theta[ik_opp]!=-theta[ik]) ik_opp++;

  return ik_opp;
}

//return the ratio between the correlation function with and without insertion
double fun_ratio(double A,double SL,double M,int t,int TH)
{
  return A+SL*(t-TH)*tanh(M*(TH-t));
}

void check_interval(int var,int min,int max)
{
  if(var>=max || var<min)
    {
      char err_mess[1024];
      sprintf(err_mess,"Asked for correlation of impossible combination, %d not in the interval: [%d,%d)\n",var,min,max);
      crash(err_mess,1);
    }
}

//read a particular two point passed as argument
void read_two_points(jack_vec *c,const char *in,int nmoms,int nmass,int im1,int im2,int ik1,int ik2,int r1,int r2,int ri)
{
  check_interval(ik1,0,nmoms);check_interval(im1,0,nmass);check_interval(r1,0,2);
  check_interval(ik2,0,nmoms);check_interval(im2,0,nmass);check_interval(r2,0,2);
  check_interval(ri,0,2);
  
  int icorr=2*(ik1+nmoms*(ik2+nmoms*(im1*2+r1+2*nmass*(im2*2+r2))))+ri;
  jack_vec_load_nazario_format(c,in,icorr);
}

//read a particular three point passed as argument
void read_three_points(jack_vec *c,const char *in,int nmoms,int nmass,int im1,int im2,int im3,int ik1,int ik2,int r1,int r2,int r3,int mu,int ri)
{
  check_interval(ik1,0,nmoms);check_interval(im1,0,nmass);check_interval(r1,0,2);
  check_interval(ik2,0,nmoms);check_interval(im2,0,nmass);check_interval(r2,0,2);
  check_interval(im3,0,nmass);check_interval(r3,0,2);
  check_interval(ri,0,2);
  check_interval(mu,0,4);

  int icorr=ik1+nmoms*(ik2+nmoms*(im1*2+r1+2*nmass*(im2*2+r2+2*nmass*(im3*2+r3))));

  FILE *file=open_file(in,"r");
  int nel=c->nel;

  if(fseek(file,icorr*2*4*sizeof(jack)*nel,SEEK_SET))
    {
      fprintf(stderr,"Error while searching for correlation %d!\n",icorr);
      exit(1);
    }

  double data_in[nel][4][2][njack+1];
  int stat=fread(data_in,sizeof(jack)*4*2*nel,1,file);
  if(stat!=1)
    {
      if(stat==EOF) crash("Error, reached EOF while reading data!\n",1);
      else
        {
          perror("Error while reading data");
          exit(1);
        }
    }
  
  for(int iel=0;iel<nel;iel++)
    for(int ijack=0;ijack<njack+1;ijack++)
      c->data[iel][ijack]=data_in[iel][mu][ri][ijack];
}

void read_P5_Vmu_P5(jack_vec *c,const char *base_path,int nmoms,int nmass,int im1,int im2,int im3,int ik1,int ik2,int r1,int r2,int r3,int mu)
{
  char path[1024];
  sprintf(path,"%s/oPVmuPo-sss_conf.1.dat",base_path);
  
  //real or imag part, according to mu
  int ri[4]={0,1,1,1};

  read_three_points(c,path,nmoms,nmass,im1,im2,im3,ik1,ik2,r1,r2,r3,mu,ri[mu]);
}

void read_P5_P5(jack_vec *c,const char *base_path,int nmoms,int nmass,int im1,int im2,int ik1,int ik2,int r1,int r2)
{
  char path[1024];
  sprintf(path,"%s/oPPo-ss_conf.1.dat",base_path);
  
  //read real part
  int ri=0;
  read_two_points(c,path,nmoms,nmass,im1,im2,ik1,ik2,r1,r2,ri);
}

void read_A0_P5(jack_vec *c,const char *base_path,int nmoms,int nmass,int im1,int im2,int ik1,int ik2,int r1,int r2)
{
  char path[1024];
  sprintf(path,"%s/oA0Po-ss_conf.1.dat",base_path);
  
  //read real part
  int ri=0;
  read_two_points(c,path,nmoms,nmass,im1,im2,ik1,ik2,r1,r2,ri);
}

void read_improved_P5_Vmu_P5(jack_vec *c,const char *base_path,int nmoms,int nmass,int im1,int im2,int im3,int ik1,int ik2,int r1,int r2,int r3,int mu,double *theta)
{
  char path[1024];
  sprintf(path,"%s/oPVmuPo-sss_conf.1.dat",base_path);
  
  jack_vec *ctemp;
  ctemp=jack_vec_malloc(c->nel);
  
  //real or imag part, according to mu
  int ri[4]={0,1,1,1};

  read_three_points(c,path,nmoms,nmass,im1,im2,im3,ik1,ik2,r1,r2,r3,mu,ri[mu]);
  read_three_points(c,path,nmoms,nmass,im1,im3,im2,opposite_theta(ik1,theta),opposite_theta(ik2,theta),r1,r2,r3,mu,ri[mu]);

  jack_vec_summassign_jack_vec(c,ctemp);
  jack_vec_prodassign_double(c,0.5);
  
  jack_vec_free(ctemp);
}

void read_improved_P5_P5(jack_vec *c,const char *base_path,int nmoms,int nmass,int im1,int im2,int ik1,int ik2,int r1,int r2,double *theta)
{
  char path[1024];
  sprintf(path,"%s/oPPo-ss_conf.1.dat",base_path);
  
  jack_vec *ctemp;
  ctemp=jack_vec_malloc(c->nel);
  
  //read real part
  int ri=0;
  read_two_points(c,path,nmoms,nmass,im1,im2,ik1,ik2,r1,r2,ri);
  read_two_points(ctemp,path,nmoms,nmass,im1,im2,opposite_theta(ik1,theta),opposite_theta(ik2,theta),r1,r2,ri);
  
  jack_vec_summassign_jack_vec(c,ctemp);
  jack_vec_prodassign_double(c,0.5);
  
  jack_vec_free(ctemp);
}

//load the mesonic three point function
void load_improved_charged_meson_three_points_P5_V0_P5(jack_vec *P5_V0_P5,const char *base_path,int nmoms,int nmass,int r,int im1,int im2,int im3,int ik1,int ik2,double *theta,int L)
{
  //load the charged "r" flavour
  int r1=r,r2=!r1,r3=r1;

  //load mu=0
  int mu=0;

  //read
  read_improved_P5_Vmu_P5(P5_V0_P5,base_path,nmoms,nmass,im1,im2,im3,ik1,ik2,r1,r2,r3,mu,theta);

  //Put the 1/spat_vol factor
  jack_vec_prodassign_double(P5_V0_P5,1.0/(L*L*L));
}

void load_improved_degenerate_charged_meson_three_points_P5_V0_P5(jack_vec *P5_V0_P5,const char *base_path,int nmoms,int nmass,int r,int im_spec,int im_valence,int ik1,int ik2,double *theta,int L)
{load_improved_charged_meson_three_points_P5_V0_P5(P5_V0_P5,base_path,nmoms,nmass,r,im_spec,im_valence,im_valence,ik1,ik2,theta,L);}

//load the mesonic two point function
void load_improved_charged_meson_two_points_P5_P5(jack_vec *P5_P5,const char *base_path,int nmoms,int nmass,int r,int im_spec,int im_valence,int ik1,double *theta,int L)
{
  //load the charged "r" flavour
  int r1=r,r2=r1;

  //load the degenerate three point
  int im1=im_spec,im2=im_valence;
  
  //load the meson
  int ik2=0;

  //read
  read_improved_P5_P5(P5_P5,base_path,nmoms,nmass,im1,im2,ik1,ik2,r1,r2,theta);
  
  //Put the -1/spat_vol factor
  jack_vec_prodassign_double(P5_P5,-1.0/(L*L*L));
}


void print_corr_to_file(const char *path,jack_vec *corr)
{
  //open the out file
  FILE *fout=open_file(path,"w");
  jack_vec_fprintf(fout,corr);  
  fclose(fout);
}

void calculate_Q2(jack Q2,jack E1,double theta1,jack E2,double theta2,int L)
{
  jack DE;
  jack_subt_jack(DE,E1,E2);
  jack_prod_jack(Q2,DE,DE);
  
  double dp=2*M_PI/L*(theta1-theta2);
  jack_subt_double(Q2,Q2,3*dp*dp);
}
