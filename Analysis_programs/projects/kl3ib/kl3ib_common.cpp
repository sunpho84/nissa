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
  if(var>max || var<min) crash("Asked for correlation of impossible combination\n",1);
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
  read_three_points(ctemp,path,nmoms,nmass,im1,im2,im3,opposite_theta(ik1,theta),opposite_theta(ik2,theta),r1,r2,r3,mu,ri[mu]);

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

