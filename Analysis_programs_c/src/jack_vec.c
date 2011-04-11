#pragma once
#include <file.c>
#include <crash.c>

typedef struct jack_vec
{
  int nel;
  jack *data;
} jack_vec;

jack_vec *jack_vec_malloc(int nel)
{
  jack_vec *i=(jack_vec*)malloc(sizeof(jack_vec));

  i->data=(jack*)malloc(sizeof(jack)*nel);
  i->nel=nel;
  
  return i;
}

void jack_vec_free(jack_vec *i)
{
  free(i->data);
  free(i);
}

void jack_vec_check_equal_nel(jack_vec *a,jack_vec *b)
{
  if(b->nel!=a->nel)
    {
      fprintf(stderr,"Number of element differs\n");
      exit(1);
    }
}

void jack_vec_load(jack_vec *v,const char *path,int i)
{
  FILE *file=open_file(path,"r");
  int nel=v->nel;

  if(fseek(file,i*sizeof(jack)*nel,SEEK_SET))
    {
      fprintf(stderr,"Error while searching for correlation %d!\n",i);
      exit(1);
    }
  
  if(fread(v->data,sizeof(jack)*nel,1,file)!=1)
    {
      fprintf(stderr,"Error while reading data!\n");
      exit(1);
    }

  fclose(file);
}

void jack_vec_load_nazario_format(jack_vec *v,const char *path,int i)
{
  FILE *file=open_file(path,"r");
  int nel=v->nel;

  if(fseek(file,2*(i/2)*sizeof(jack)*nel,SEEK_SET))
    {
      fprintf(stderr,"Error while searching for correlation %d!\n",i);
      exit(1);
    }

  double data_in[nel][2][njack+1];
  int stat=fread(data_in,sizeof(jack)*2*nel,1,file);
  if(stat!=1)
    {
      if(stat==EOF) crash("Error, reached EOF while reading data!\n",1);
      else
	{
	  perror("Error while reading data");
	  exit(1);
	}
    }

  int ri=i%2;
  for(int i=0;i<nel;i++)
    for(int ijack=0;ijack<njack+1;ijack++)
      v->data[i][ijack]=data_in[i][ri][ijack];
  
  fclose(file);
}

void jack_vec_fprintf(FILE *file,jack_vec *v)
{
  int nel=v->nel;
  
  for(int iel=0;iel<nel;iel++) fprintf(file,"%d %g %g\n",iel,v->data[iel][njack],jack_error(v->data[iel]));
}

void jack_vec_print_to_file(const char *path,jack_vec *corr)
{
  //open the out file
  FILE *fout=open_file(path,"w");
  jack_vec_fprintf(fout,corr);  
  fclose(fout);
}

void jack_vec_summ_jack_vec(jack_vec *c,jack_vec *a,jack_vec *b)
{
  int nel=a->nel;
  jack_vec_check_equal_nel(a,b);
  jack_vec_check_equal_nel(a,c);
  
  for(int iel=0;iel<nel;iel++)
    for(int ijack=0;ijack<njack+1;ijack++)
      c->data[iel][ijack]=a->data[iel][ijack]+b->data[iel][ijack];
}

void jack_vec_summassign_jack_vec(jack_vec *c,jack_vec *b)
{jack_vec_summ_jack_vec(c,c,b);}

void jack_vec_subt_jack_vec(jack_vec *c,jack_vec *a,jack_vec *b)
{
  int nel=a->nel;
  jack_vec_check_equal_nel(a,b);
  jack_vec_check_equal_nel(a,c);
  
  for(int iel=0;iel<nel;iel++)
    for(int ijack=0;ijack<njack+1;ijack++)
      c->data[iel][ijack]=a->data[iel][ijack]-b->data[iel][ijack];
}

void jack_vec_subtassign_jack_vec(jack_vec *c,jack_vec *b)
{jack_vec_subt_jack_vec(c,c,b);}

void jack_vec_prod_double(jack_vec *c,jack_vec *a,double b)
{
  int nel=a->nel;
  jack_vec_check_equal_nel(a,c);

  for(int iel=0;iel<nel;iel++)
    for(int ijack=0;ijack<njack+1;ijack++)
      c->data[iel][ijack]=a->data[iel][ijack]*b;
}

void jack_vec_prodassign_double(jack_vec *a,double b)
{jack_vec_prod_double(a,a,b);}

void jack_vec_fracassign_double(jack_vec *a,double b)
{jack_vec_prod_double(a,a,1/b);}

void jack_vec_prod_jack(jack_vec *c,jack_vec *a,jack b)
{
  int nel=a->nel;
  jack_vec_check_equal_nel(a,c);
  
  for(int iel=0;iel<nel;iel++)
    for(int ijack=0;ijack<njack+1;ijack++)
      c->data[iel][ijack]=a->data[iel][ijack]*b[ijack];
}

void jack_vec_prodassign_jack(jack_vec *c,jack b)
{jack_vec_prod_jack(c,c,b);}

void jack_vec_prod_jack_vec(jack_vec *c,jack_vec *a,jack_vec *b)
{
  int nel=a->nel;
  jack_vec_check_equal_nel(a,b);
  jack_vec_check_equal_nel(a,c);

  for(int iel=0;iel<nel;iel++)
    for(int ijack=0;ijack<njack+1;ijack++)
      c->data[iel][ijack]=a->data[iel][ijack]*b->data[iel][ijack];
}

void jack_vec_prodassign_jack_vec(jack_vec *c,jack_vec *b)
{jack_vec_prod_jack_vec(c,c,b);}

void jack_vec_average(jack_vec *corr_av,jack_vec **corr,int n)
{
  int nel=corr_av->nel;
  for(int i=0;i<n;i++)
    if(corr[i]->nel!=nel)
      {
	fprintf(stderr,"Correlation %d to average differs from other\n",i);
	exit(1);
      }
  
  memcpy(corr_av->data,corr[0]->data,sizeof(jack)*nel);
  
  for(int i=1;i<n;i++) jack_vec_summassign_jack_vec(corr_av,corr[i]);
  jack_vec_prodassign_double(corr_av,1.0/n);
}

void jack_vec_simmetrize(jack_vec *out,jack_vec *in,int si)
{
  int T=in->nel,TH=out->nel;
  jack *outd=out->data,*ind=in->data;
  
  if(T%2!=0)
    {
      fprintf(stderr,"Error in simmetrization, input of odd elements\n");
      exit(1);
    }

  if(TH!=T/2)
    {
      fprintf(stderr,"Error in simmetrization, number of elemenets of output: %d, of input: %d\n",TH,T);
      exit(1);
    }

  for(int t1=0;t1<TH;t1++)
    {
      int t2=T-t1;
      if(t2==T) t2=0;
      
      for(int ij=0;ij<njack+1;ij++)
	if(si==1) outd[t1][ij]=(ind[t1][ij]+ind[t2][ij])/2;
	else      outd[t1][ij]=(ind[t1][ij]-ind[t2][ij])/2;
    }  
}

void jack_vec_frac_jack_vec(jack_vec *c,jack_vec *a,jack_vec *b)
{
  int nel=a->nel;
  jack_vec_check_equal_nel(a,b);
  jack_vec_check_equal_nel(a,c);

  for(int iel=0;iel<nel;iel++)
    for(int ijack=0;ijack<njack+1;ijack++)
      c->data[iel][ijack]=a->data[iel][ijack]/b->data[iel][ijack];
}

void jack_vec_fracassign_jack_vec(jack_vec *c,jack_vec *b)
{jack_vec_frac_jack_vec(c,c,b);}

void jack_vec_sqrt_jack_vec(jack_vec *c,jack_vec *a)
{
  int nel=a->nel;
  jack_vec_check_equal_nel(a,c);

  for(int iel=0;iel<nel;iel++)
    for(int ijack=0;ijack<njack+1;ijack++)
      c->data[iel][ijack]=sqrt(a->data[iel][ijack]);
}

void jack_vec_frac_jack(jack_vec *c,jack_vec *a,jack b)
{
  int nel=a->nel;
  jack_vec_check_equal_nel(a,c);

  for(int iel=0;iel<nel;iel++)
    for(int ijack=0;ijack<njack+1;ijack++)
      c->data[iel][ijack]=a->data[iel][ijack]/b[ijack];
}

void jack_vec_fracassign_jack(jack_vec *c,jack b)
{jack_vec_frac_jack(c,c,b);}

void jack_frac_jack_vec(jack_vec *c,jack a,jack_vec *b)
{
  int nel=b->nel;
  jack_vec_check_equal_nel(b,c);

  for(int iel=0;iel<nel;iel++)
    for(int ijack=0;ijack<njack+1;ijack++)
      c->data[iel][ijack]=a[ijack]/b->data[iel][ijack];
}

