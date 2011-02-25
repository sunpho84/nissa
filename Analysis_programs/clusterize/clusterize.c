#include <stdio.h>
#include <string.h>
#include <stdlib.h>

typedef char string[1024];

int nlevel;
string *var_name;
int *var_value_tot;
int *required_header;
string **var_values;
string *level_header_template;

FILE *open_file(const string path,const string mod)
{
  FILE *out=fopen(path,mod);

  if(out==NULL)
    {
      fprintf(stderr,"Error opening file %s for %s\n",path,mod);
      exit(1);
    }
  
  return out;
}

void expect_string(FILE *fin,const string exp)
{
  string read;
  int nscan=fscanf(fin,"%s",read);
  if(nscan!=1||strcmp(exp,read)!=0)
    {
      if(nscan==EOF) fprintf(stderr,"Error,reached file end while waiting for %s\n",exp);
      else           fprintf(stderr,"Error, read %s instead than %s\n",read,exp);
      exit(1);
    }
}

void read_file(char *out,FILE *fin,const string what,const string varname)
{
  int nscan=fscanf(fin,what,out);
  if(nscan!=1)
    {
      if(nscan==EOF) fprintf(stderr,"Error,reached file end while reading %s\n",varname);
      else           fprintf(stderr,"Error, not enough data while reading %s\n",varname);
      exit(1);
    }
}

void read_file_expecting(char *out,FILE *fin,const string what,const string varname)
{
  expect_string(fin,varname);
  read_file(out,fin,what,varname);
}

void prepare_header(string header,string template,int *level_entry)
{
  memcpy(header,template,1024);

  for(int ilevel=0;ilevel<nlevel;ilevel++)
    {
      int lheader=strlen(header);
      
      //prepare the tag
      string tag="[";
      strcat(tag,var_name[ilevel]);
      strcat(tag,"]");
      int tag_length=strlen(tag);
      
      //find the position of the tag
      char *tag_start_pos=strstr(header,tag);
      int lfirst_part=(int)(tag_start_pos-header);
      char *tag_end_pos=tag_start_pos+tag_length;

      if(lfirst_part>=0 && lfirst_part<lheader)
	{
	  string temp_header="";
	  strncat(temp_header,header,lfirst_part);
	  strcat(temp_header,var_values[ilevel][level_entry[ilevel]]);
	  strcat(temp_header,tag_end_pos);
	  memcpy(header,temp_header,sizeof(string));
	}
    }
  
  int lheader=strlen(header);
  for(int i=0;i<lheader;i++) if(header[i]=='\'') header[i]=' ';
}

int main(int narg,char **arg)
{
  if(narg<2)
    {
      fprintf(stderr,"Error, use: %s run_input\n",arg[0]);
      exit(1);
    }

  FILE *frun_input=open_file(arg[1],"r");
  
  //read the analysis template
  string analysis_template_path;
  read_file_expecting(analysis_template_path,frun_input,"%s","analysis_template_path");
  FILE *fatempl=open_file(analysis_template_path,"r");
  read_file_expecting((char*)&nlevel,fatempl,"%d","nlevel");
  //allocate room for the header format
  var_name=(string*)malloc(nlevel*sizeof(string));
  var_value_tot=(int*)malloc(nlevel*sizeof(int));
  var_values=(string**)malloc(nlevel*sizeof(string*));
  level_header_template=(string*)malloc(nlevel*sizeof(string));
  required_header=(int*)malloc(nlevel*sizeof(int));
  //read the name of the looping variables and their values
  for(int ilevel=0;ilevel<nlevel;ilevel++)
    {
      read_file_expecting(var_name[ilevel],fatempl,"%s","var_name");
      expect_string(frun_input,var_name[ilevel]);
      read_file_expecting((char*)&(var_value_tot[ilevel]),frun_input,"%d","var_value_tot");

      var_values[ilevel]=(string*)malloc(var_value_tot[ilevel]*sizeof(string));
      expect_string(frun_input,"list_of_var_values");
      for(int ival=0;ival<var_value_tot[ilevel];ival++)
	read_file(var_values[ilevel][ival],frun_input,"%s","var_value");
      
      read_file_expecting(level_header_template[ilevel],fatempl,"%s","header");
      required_header[ilevel]=strcmp(level_header_template[ilevel],"NO");
    }
  fclose(fatempl);
  //read T
  int T;
  read_file_expecting((char*)&T,frun_input,"%d","T");
  //read the path template
  string path_template;
  read_file_expecting(path_template,frun_input,"%s","path_template");
  //read the number of files to open
  int nfiles;
  read_file_expecting((char*)&nfiles,frun_input,"%d","nfiles");
  //open all the input files, reading each entry
  expect_string(frun_input,"list_of_files");
  FILE *fdata[nfiles];
  for(int ifile=0;ifile<nfiles;ifile++)
    {
      string chunk;
      read_file(chunk,frun_input,"%s","particular file name");
      string file_path;
      sprintf(file_path,path_template,chunk);
      fdata[ifile]=open_file(file_path,"r");
    }
  //read the output file
  string outpath;
  read_file_expecting(outpath,frun_input,"%s","outfile");
  FILE *fout=open_file(outpath,"w");

  /////////////////////////////////////////////////////////////////////////////////

  //calculate the number of combo on which to loop
  int ncombo=1;
  for(int ilevel=0;ilevel<nlevel;ilevel++) ncombo*=var_value_tot[ilevel];
  
  //loop over all the combo
  for(int icombo=0;icombo<ncombo;icombo++)
    {
      int level_entry[nlevel];
      memset(level_entry,0,sizeof(int)*nlevel);
      
      int temp_icombo=icombo;
      for(int ilevel=nlevel-1;ilevel>=0;ilevel--)
	{
	  level_entry[ilevel]=temp_icombo%var_value_tot[ilevel];
	  temp_icombo/=var_value_tot[ilevel];
	}
      
      //now read all the required header
      for(int ilevel=0;ilevel<nlevel;ilevel++)
	if(required_header[ilevel] && (ilevel==nlevel-1||level_entry[ilevel+1]==0))
	  {
	    string header;
	    prepare_header(header,level_header_template[ilevel],level_entry);
	    
	    for(int ifile=0;ifile<nfiles;ifile++)
	      {
		string read_header;
		char *read=fgets(read_header,1024,fdata[ifile]);
		if(read!=read_header)
		  {
		    fprintf(stderr,"Error while reading header: %s\n",header);
		    exit(1);
		  }
		if(strcmp(read,read_header))
		  {
		    fprintf(stderr,"Error while reading header: %s, obtained: %s\n",header,read_header);
		    exit(1);
		  }
	      }
	  }
      
      //now read a whole correlation
      double data[nfiles][T][2];
      for(int ifile=0;ifile<nfiles;ifile++)
	for(int t=0;t<T;t++)
	  for(int ri=0;ri<2;ri++)
	    read_file((char*)&data[ifile][t][ri],fdata[ifile],"%g","data");
    }

  return 0;
}
