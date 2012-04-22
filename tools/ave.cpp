#include <stdio.h>
#include <stdlib.h>

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
    if(narg<4)
      {	
	fprintf(stderr,"Error,use %s n listin out\n",arg[0]);
	exit(1);
      }
    int nin=atoi(arg[1]);
    if(narg!=3+nin)
      {
        fprintf(stderr,"Error, expected %d in listin,only %d present\n",nin,narg-3);
        exit(1);
      }
    
    FILE *f[nin];
    for(int i=0;i<nin;i++)
    {
	f[i]=fopen(arg[i+2],"r");
	if(f[i]==NULL)
	{
	    fprintf(stderr,"Error opening file %s\n",arg[i+2]);
	    exit(1);
	}
    }
    
    FILE *out=fopen(arg[narg-1],"w");
    if(out==NULL)
    {
	fprintf(stderr,"Error opening file %s\n",arg[narg-1]);
	exit(1);
    }

    
    double er;
    do
    {
	double a=0,t;
	er=0;
	for(int i=0;i<nin;i++)
	{
	    er+=fread(&t,sizeof(double),1,f[i]);
	    //change(&t,&t,1);
	    a+=t;
	}
	a/=nin;
	//printf("%lg\n",a);
	//change(&a,&a,1);
	if(er==nin) fwrite(&a,sizeof(double),1,out);
    }
    while(er==nin);
    
    return 0;
}
	
