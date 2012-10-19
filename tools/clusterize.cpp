#include <stdio.h>
#include <string.h>
#include <stdlib.h>

int check_endianess()
{
  int big_endian=1;
  big_endian=(int)(*(char*)(&big_endian));
  return big_endian;
}

int count_corr(char *path,const char *line_to_find)
{
  int ncorr=0;
  
  FILE *fin=fopen(path,"r");
  if(fin==NULL)
    {
      fprintf(stderr,"Error opening file: %s\n",path);
      exit(1);
    }

    char line[1024];
    while(fgets(line,1024,fin)==line)
	if(strcmp(line,line_to_find)==0)
	  ncorr++;
    
    fclose(fin);
    
    printf("found %d entries\n",ncorr);
    
    return ncorr;
}

void change(double *dest,double *sour,int ndoubles)
{
  char *cdest,*csour;
  char temp;
  
  for(int idouble=0;idouble<ndoubles;idouble++)
    {
      cdest=(char*)(dest+idouble);
      csour=(char*)(sour+idouble);
      
      temp=csour[7];
      csour[7]=cdest[0];
      cdest[0]=temp;
      
      temp=csour[6];
      csour[6]=cdest[1];
      cdest[1]=temp;
      
      temp=csour[5];
      csour[5]=cdest[2];
      cdest[2]=temp;
      
      temp=csour[4];
      csour[4]=cdest[3];
      cdest[3]=temp;
    }
}

int main(int narg,char **arg)
{
  int end=check_endianess();
  
  if(narg<9)
    {
      fprintf(stderr,"Use: %s base_in base_gauge ngauge njack out L corr mode[0=reim,1=re]\n",arg[0]);
      exit(0);
    }
  
  int base_gauge=atoi(arg[2]);
  int nconf_teo=atoi(arg[3]);
  int njack=atoi(arg[4]);
  int clust_size=nconf_teo/njack;
  int nconf=clust_size*njack;
  int T=atoi(arg[6])*2;
  int mode=atoi(arg[8]);
  
  char line_to_find[20];
  sprintf(line_to_find," # %s\n",arg[7]);
  
  char path[1024];
  sprintf(path,arg[1],base_gauge);
  int ncorr=count_corr(path,line_to_find);
  
  //load data
  double ****data=(double****)malloc(ncorr*sizeof(double***));
  for(int icorr=0;icorr<ncorr;icorr++)
    {
      data[icorr]=(double***)malloc(2*sizeof(double**));
      for(int ri=0;ri<2;ri++)
	{
	  data[icorr][ri]=(double**)malloc(T*sizeof(double*));
	  for(int t=0;t<T;t++) data[icorr][ri][t]=(double*)malloc(nconf*sizeof(double));
	}
      }
  
  int targ_conf=0;
  for(int iconf=0;iconf<nconf;iconf++)
    {
      FILE *fin;
      
      printf("%d\n",iconf);
	
      //open file
      do
	{
	  sprintf(path,arg[1],targ_conf+base_gauge);
	  fin=fopen(path,"r");
	  if(fin==NULL)
	    {
	      fprintf(stderr,"Error, couldn't open file: %s, skipping\n",path);
	      targ_conf++;
	      if(targ_conf>=nconf_teo)
		{
		  fprintf(stderr,"Error, finished the available confs\n");
		  exit(1);
		}
	    }
	  
	}
      while(fin==NULL);
      targ_conf++;
      
      //search the corr
      for(int icorr=0;icorr<ncorr;icorr++)
	{
	  char line[1024],*out;
	  do out=fgets(line,1024,fin);
	  while(out==line && (strcmp(line,line_to_find)!=0));
	  
	  //read
	  for(int t=0;t<T;t++)
	    {
	      int r=fscanf(fin,"%lg %lg",&(data[icorr][0][t][iconf]),&(data[icorr][1][t][iconf]));
	      if(r!=2)
		{
		  fprintf(stderr,"Error reading from file %d time: %d\n",iconf,t);
		  exit(1);
		}
	    }
	}
      
      //close file
      fclose(fin);
    }
  
  FILE *fout=fopen(arg[5],"w");
  if(fout==NULL)
    {
      fprintf(stderr,"Error, couldn't open for output file: %s\n",arg[5]);
      exit(1);
    }
  
  //clusterize
  for(int icorr=0;icorr<ncorr;icorr++)
    {
      double jack[2][T][njack+1];
      for(int ri=0;ri<2;ri++)
	for(int t=0;t<T;t++)
	  {
	    double clus[njack+1];
	    memset(clus,0,sizeof(double)*(njack+1));
	    
	    for(int iconf=0;iconf<nconf;iconf++) clus[iconf/clust_size]+=data[icorr][ri][t][iconf];
	    for(int ijack=0;ijack<njack;ijack++) clus[njack]+=clus[ijack];
	    for(int ijack=0;ijack<njack;ijack++) jack[ri][t][ijack]=(clus[njack]-clus[ijack])/(nconf-nconf/njack);
	    jack[ri][t][njack]=clus[njack]/nconf;
	  }
      
      //write
      if(end==0) change((double*)jack,(double*)jack,T*2*(njack+1));
      if(mode==0) fwrite(jack,sizeof(double),T*2*(njack+1),fout);
      else        fwrite(jack,sizeof(double),T*(njack+1),fout);
    }
  
  fclose(fout);
  
  return 0;
}
