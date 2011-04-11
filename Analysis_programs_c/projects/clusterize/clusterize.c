#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include <analysis_include.h>

typedef char string[1024];

int nlevel;
string *var_name;
int *var_value_tot;
int *required_header;
string **var_values;
string *level_header_templ;

void prepare_header(string header,string templ,int *level_entry)
{
  memcpy(header,templ,1024);

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
  strcat(header,"\n");
}

int main(int narg,char **arg)
{
  if(narg<2)
    {
      fprintf(stderr,"Error, use: %s run_input\n",arg[0]);
      exit(1);
    }
  
  int big_endian=check_endianess();
  
  FILE *frun_input=open_file(arg[1],"r");
  
  //read the analysis templ
  string analysis_templ_path;
  read_formatted_from_file_expecting(analysis_templ_path,frun_input,"%s","analysis_templ_path");
  FILE *fatempl=open_file(analysis_templ_path,"r");
  read_formatted_from_file_expecting((char*)&nlevel,fatempl,"%d","nlevel");
  //allocate room for the header format
  var_name=(string*)malloc(nlevel*sizeof(string));
  var_value_tot=(int*)malloc(nlevel*sizeof(int));
  var_values=(string**)malloc(nlevel*sizeof(string*));
  level_header_templ=(string*)malloc(nlevel*sizeof(string));
  required_header=(int*)malloc(nlevel*sizeof(int));
  //read the name of the looping variables and their values
  for(int ilevel=0;ilevel<nlevel;ilevel++)
    {
      read_formatted_from_file_expecting(var_name[ilevel],fatempl,"%s","var_name");
      expect_string_from_file(frun_input,var_name[ilevel]);
      read_formatted_from_file_expecting((char*)&(var_value_tot[ilevel]),frun_input,"%d","var_value_tot");

      var_values[ilevel]=(string*)malloc(var_value_tot[ilevel]*sizeof(string));
      expect_string_from_file(frun_input,"list_of_var_values");
      for(int ival=0;ival<var_value_tot[ilevel];ival++)
	read_formatted_from_file(var_values[ilevel][ival],frun_input,"%s","var_value");
      
      read_formatted_from_file_expecting(level_header_templ[ilevel],fatempl,"%s","header");
      required_header[ilevel]=strcmp(level_header_templ[ilevel],"NO");
    }
  fclose(fatempl);
  //read T
  int T;
  read_formatted_from_file_expecting((char*)&T,frun_input,"%d","T");
  //read the output file
  string outpath;
  read_formatted_from_file_expecting(outpath,frun_input,"%s","outfile");
  FILE *fout=open_file(outpath,"w");
  //read the number of jacknives
  int njack_read;
  read_formatted_from_file_expecting((char*)&njack_read,frun_input,"%d","njack");
  if(njack!=njack_read)
    {
      fprintf(stderr,"Error, number of jacknives fixed to: %d, read: %d.\n",njack,njack_read);
      exit(1);
    }
  //read the path templ
  string path_templ;
  read_formatted_from_file_expecting(path_templ,frun_input,"%s","path_templ");
  //read the number of files to open
  int tfiles;
  read_formatted_from_file_expecting((char*)&tfiles,frun_input,"%d","nfiles");
  //calculate cluster size
  int clus_size=tfiles/njack;
  int nfiles=clus_size*njack;
  printf("Njack: %d, cluster size: %d\n",njack,clus_size);
  //open all the input files, reading each entry
  expect_string_from_file(frun_input,"list_of_files");
  FILE *fdata[nfiles];
  for(int ifile=0;ifile<nfiles;ifile++)
    {
      string chunk;
      read_formatted_from_file(chunk,frun_input,"%s","particular file name");
      string file_path;
      sprintf(file_path,path_templ,chunk);
      fdata[ifile]=open_file(file_path,"r");
      //printf("file %d = %s\n",ifile,file_path);
    }

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
	    prepare_header(header,level_header_templ[ilevel],level_entry);
	    
	    for(int ifile=0;ifile<nfiles;ifile++)
	      {
		string read_header;
		int nchar;
		do
		  {
		    char *read=fgets(read_header,1024,fdata[ifile]);
		    if(read!=read_header)
		      {
			fprintf(stderr,"Error while reading header: %s\n",header);
			exit(1);
		      }
		    
		    nchar=0;
		    do if(*read!=' ' && *read!='\n' && *read!='\0') nchar++;
		    while(*(read++)!='\0' && nchar==0);
		  }
		while(nchar==0);
		  
		  if(strcmp(read_header,header))
		    {
		      fprintf(stderr,"Error while reading header: '%s', obtained: '%s'\n",header,read_header);
		      exit(1);
		    }
	      }
	  }
      
      //now read a whole correlation
      double data[nfiles][T][2];
      for(int ifile=0;ifile<nfiles;ifile++)
	for(int t=0;t<T;t++)
	  for(int ri=0;ri<2;ri++)
	    read_formatted_from_file((char*)&data[ifile][t][ri],fdata[ifile],"%lg","data");
      
      //clusterize
      double clus[2][T][njack+1];
      memset(clus,0,sizeof(double)*2*T*(njack+1));
      for(int ri=0;ri<2;ri++)
	for(int t=0;t<T;t++)
	  {
	    for(int ifile=0;ifile<nfiles;ifile++)
	      {
		int iclus=ifile%njack;
		clus[ri][t][iclus]+=data[ifile][t][ri];
	      }
	    for(int iclus=0;iclus<njack;iclus++) clus[ri][t][njack]+=clus[ri][t][iclus];
	    for(int ijack=0;ijack<njack;ijack++) clus[ri][t][ijack]=(clus[ri][t][njack]-clus[ri][t][ijack])/(nfiles-clus_size);
	    clus[ri][t][njack]/=nfiles;
	  }
      
      if(big_endian==0) doubles_to_doubles_changing_endianess((double*)clus,(double*)clus,2*T*(njack+1));
      
      for(int t=0;t<T;t++) printf("%d %lg %lg\n",t,clus[0][t][njack],clus[1][t][njack]);
      
      if(fwrite(clus,sizeof(double),2*T*(njack+1),fout)!=2*T*(njack+1))
	{
	  fprintf(stderr,"Error while writing double!\n");
	  exit(1);
	}

      if(icombo%100==0) printf("Completed correlation: %d\n",icombo);
    }

  return 0;
}
