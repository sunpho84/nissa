int little_endian;
FILE *input_global;

//check the endianess of the machine
void check_endianess()
{
  little_endian=1;
  little_endian=(int)(*(char*)(&little_endian));
}

//fucking tool to revert the endianess of doubles
void doubles_to_doubles_changing_endianess(double *dest,double *sour,int ndoubles)
{
  char *cdest,*csour;
  char temp;
  
  if(dest==sour)
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
  else
    for(int idouble=0;idouble<ndoubles;idouble++)
      {
        cdest=(char*)(dest+idouble);
        csour=(char*)(sour+idouble);
        
        cdest[0]=csour[7];
        cdest[1]=csour[6];
        cdest[2]=csour[5];
        cdest[3]=csour[4];
        cdest[4]=csour[3];
        cdest[5]=csour[2];
        cdest[6]=csour[1];
        cdest[7]=csour[0];

      }
}

void CRASH(const char *templ,...)
{
  va_list ap;
  va_start(ap,templ);
  
  char mess[1024];
  vsprintf(mess,templ,ap);
  va_end(ap);
  
  fprintf(stderr,"%s\n",mess);  
  
  exit(1);
}

//read a token from file
int read_next_token(char *tok)
{
  int ok=fscanf(input_global,"%s",tok);

  return ok;
}

//check whether a token starts or end a comment
int check_tok_starts_comment(char *tok)
{return strncasecmp(tok,"/*",2)==0;}
int check_tok_ends_comment(char *tok)
{return strncasecmp(tok+strlen(tok)-2,"*/",2)==0;}

//read up to "*/"
void read_up_to_end_of_comment()
{
  char tok[1024];
  
  do
    {
      int ok=read_next_token(tok);
      if(ok!=1) CRASH("reached end of file without finding end of comment");
    }
  while(check_tok_ends_comment(tok)==0);
}

int read_var_catcherr(char *out,const char *par,int size_of)
{
  char tok[1024];
  int ok,comment_found;
  
  do
    {
      //read the token
      ok=read_next_token(tok);
      
      //check if token comment
      if(check_tok_starts_comment(tok))
        {
          comment_found=1;
          //if the comments do not ends itself, read up to finding its end
          if(!check_tok_ends_comment(tok)) read_up_to_end_of_comment();
        }
      else comment_found=0;
    }
  while(comment_found);
  
  //parse the token
  if(ok==1)
    ok=(sscanf(tok,par,out)>0);
  
  return ok;
}

void read_var(char *out,const char *par,int size_of)
{if(!read_var_catcherr(out,par,size_of)) CRASH("Couldn't read from input file!!!");}

//Read an integer from the file
void read_int(int *out)
{read_var((char*)out,"%d",sizeof(int));}

//Read a double from the file
void read_double(double *out)
{read_var((char*)out,"%lg",sizeof(double));}

//Read a string from the file
void read_str(char *str,int length)
{read_var(str,"%s",length);}

//Read a string from the file and check against the argument
void expect_str(const char *exp_str)
{
  char obt_str[1024];
  
  read_str(obt_str,1024);
  
  if(strcasecmp(exp_str,obt_str)!=0) CRASH("Error, expexcted '%s' in input file, obtained: '%s'",exp_str,obt_str);
}

//Read an int checking the tag
void read_str_int(const char *exp_str,int *in)
{
  expect_str(exp_str);
  read_int(in);
}

//Read a double checking the tag
void read_str_double(const char *exp_str,double *in)
{
  expect_str(exp_str);
  read_double(in);
}

//Read a string checking the tag
void read_str_str(const char *exp_str,char *in,int length)
{
  expect_str(exp_str);
  read_str(in,length);
}

FILE* open_file(const char *outfile,const char *mode)
{
  FILE *fout=NULL;
  
  fout=fopen(outfile,mode);
  if(fout==NULL) CRASH("Couldn't open the file: %s for mode: %s",outfile,mode);
  
  return fout;
}

//check if a file exists
int file_exists(const char *path)
{
  int status=1;
  
  FILE *f=fopen(path,"r");
  if(f!=NULL)
    {
      status=1;
      fclose(f);
    }
  else status=0;

  return status;
}

void open_input(char *input_path)
{
  input_global=fopen(input_path,"r");
  if(input_global==NULL) CRASH("File '%s' not found",input_path);
}
