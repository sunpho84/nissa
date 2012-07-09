#include "nissa.h"

int read_size;
double *data=NULL,*ave_data=NULL;

//read the file
void read_file(char *inpath)
{
  //open file
  FILE *fin=open_file(inpath,"r");
  
  //read
  int n=fread((void*)data,sizeof(double),read_size,fin);
  if(n!=read_size) crash("while reading %s received %d",inpath,n);
  
  //change endianess if needed
  if(!little_endian)
    doubles_to_doubles_changing_endianess(data,data,read_size);
  
  //summ
  for(int i=0;i<read_size;i++) ave_data[i]+=data[i];
}

//write the file
void write_file(char *outpath)
{
  //open file
  FILE *fout=open_file(outpath,"w");
  
  //change endianess if needed
  if(!little_endian)
    doubles_to_doubles_changing_endianess(ave_data,ave_data,read_size);
  
  //read
  int n=fwrite((void*)ave_data,sizeof(double),read_size,fout);
  if(n!=read_size) crash("while writing %s received %d",outpath,n);
}

//get next token
int read_next(char *tok,char *&line)
{
  int n=sscanf(line,"%s",tok);
  if(n>0)
    {
      //printf("read %d chars\n",n);
      int l=strlen(tok);
      line=strstr(line,tok);
      if(line!=NULL) line+=l;
      else crash("got null!");
    }
  else line=NULL;
  
  return n;
}

//read the number of double in a file
int get_file_ndouble(char *path)
{
  //open the file and go to the end
  FILE *fin=open_file(path,"r");
  fseek(fin,0,SEEK_END);
  
  //get position and close
  int end=ftell(fin);
  fclose(fin);
  
  //convert to ndouble
  if(end%sizeof(double)) crash("file %s contains wrong data type (size=%d, not multiple of %d)",path,end,sizeof(double));
  end/=sizeof(double);
  
  return end;
}

//average all inputs and write to output
void average(char *read_line)
{
  //get outpath
  char outpath[1024];
  if(read_next(outpath,read_line)<0) crash("reading output path");
  
  //allocated size
  int allo_size;
  
  //one by one get inpath
  int nin=0;
  do
    {
      char inpath[1024];
      if(read_next(inpath,read_line)>0)
	{
	  //get file size
	  read_size=get_file_ndouble(inpath);
	  
	  //if first file
	  if(nin==0)
	    {
	      allo_size=read_size;
	      //allocate data
	      data=nissa_malloc("data",allo_size,double);
	      ave_data=nissa_malloc("ave_data",allo_size,double);
	      memset(ave_data,0,sizeof(double)*allo_size);
	    }
	  else
	    //check agreement of all files
	    if(read_size!=allo_size) crash("file %s contains %d double, previous ones contains %d",inpath,read_size,allo_size);
	  
	  //read the file and summ it
	  read_file(inpath);
	}
      else if(nin==0) crash("while reading first inpath");
      nin++;
    }
  while(read_line!=NULL);
  
  //average
  for(int i=0;i<read_size;i++) ave_data[i]/=nin;
  
  //write file
  write_file(outpath);
  
  nissa_free(data);
  nissa_free(ave_data);
}

int main(int narg,char **arg)
{
  init_nissa();
  
  if(narg<2) crash("use %s file",arg[0]);

  //Open input
  FILE *finput;
  if(arg[1][0]=='-') finput=stdin;
  else finput=fopen(arg[1],"r");
  
  //read one by one every line
  char read_line[1024];
  while(fgets(read_line,1024,finput))
    {
      //parse the line
      int len=strlen(read_line);
      if(len>1) average(read_line);
    }
  
  close_nissa();
  
  return 0;
}
