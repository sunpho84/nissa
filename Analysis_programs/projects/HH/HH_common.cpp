#include "include.h"
#include "../../nf2/common_pars.cpp"

void read_ensemble_pars(char *base_path,int &T,int &ibeta,int &nmass,double *&mass,double &sea_mass,const char *data_list_file)
{
  FILE *input=open_file(data_list_file,"r");
  
  read_formatted_from_file_expecting(base_path,input,"%s","base_path");
  read_formatted_from_file_expecting((char*)&T,input,"%d","T");
  read_formatted_from_file_expecting((char*)&ibeta,input,"%d","njack"); //dum
  read_formatted_from_file_expecting((char*)&ibeta,input,"%d","Beta");
  
  read_formatted_from_file_expecting((char*)&nmass,input,"%d","nmass");

  expect_string_from_file(input,"mass_list");
  mass=(double*)malloc(sizeof(double)*nmass);
  for(int imass=0;imass<nmass;imass++) read_formatted_from_file((char*)&(mass[imass]),input,"%lg","mass");
  read_formatted_from_file_expecting((char*)&sea_mass,input,"%lg","sea_mass");
  
  fclose(input);
}

int icombo(int im1,int im2,int nmass)
{
 return im1*nmass+im2;
}

//load all data
void load_ensembles_list(char **&base_corrs_path,char **&ens_name,int &nens,int *&T,int *&ibeta,int *&nmass,double **&mass,double *&sea_mass,const char *ens_list_path)
{
  FILE *ens_list_file=open_file(ens_list_path,"r");
  
  //load the number of ensembles
  read_formatted_from_file_expecting((char*)(&nens),ens_list_file,"%d","nens");
  
  //base address of all correlations
  char base_ens_path[1024];
  read_formatted_from_file_expecting((char*)(&base_ens_path),ens_list_file,"%s","base_ens_path");
  
  //ensemble parameters
  base_corrs_path=(char**)malloc(nens*sizeof(char*));
  ens_name=(char**)malloc(nens*sizeof(char*));
  T=(int*)malloc(nens*sizeof(int));
  ibeta=(int*)malloc(nens*sizeof(int));
  nmass=(int*)malloc(nens*sizeof(int));
  mass=(double**)malloc(nens*sizeof(double*));
  sea_mass=(double*)malloc(nens*sizeof(double));
  
  //Loop over all ensemble. Ensemble specification is supposed to be stored in a file named [base_ens_path]/[ens_path]/data_list
  for(int iens=0;iens<nens;iens++)
    {
      //allocate ensemble name and correlations path
      base_corrs_path[iens]=(char*)malloc(1024);
      ens_name[iens]=(char*)malloc(1024);
      
      //load ensemble name and parameters
      read_formatted_from_file(ens_name[iens],ens_list_file,"%s","ens_name");
      read_ensemble_pars(base_corrs_path[iens],T[iens],ibeta[iens],nmass[iens],mass[iens],sea_mass[iens],combine("%s/%s/data_list",base_ens_path,ens_name[iens]).c_str());
    }
  
  fclose(ens_list_file);
}

void boot_from_jack(boot &out,jack in,int *iboot_jack)
{
  int nboot=out.nboot;
  int njack=in.njack;
  for(int iboot=0;iboot<nboot;iboot++) out.data[iboot]=in.data[iboot_jack[iboot]];
  out.data[nboot]=in.data[njack];
}

//load all data
void load_all_ensembles_M(bvec *&M,int &nens,int *&T,int *&ibeta,int *&nmass,const char *base_M_path,const char *obs_name,char **ens_name,char **base_corrs_path,const char *M_name)
{
  init_latpars();
  
  //allocate room for m
  M=(bvec*)malloc(nens*sizeof(bvec));
  
  //Loop over ensembles. Data is supposed to be stored in a file named [base_M_path]/[obs_name]/[ens_name]/results
  for(int iens=0;iens<nens;iens++)
    {
      char M_path[1024];
      sprintf(M_path,"%s/%s/%s/%s",base_M_path,obs_name,ens_name[iens],M_name);
      
      cout<<"Reading ensemble: "<<ens_name[iens]<<" from file: "<<M_path<<endl;
      
      int ncombo=nmass[iens]*nmass[iens];
      M[iens].create(ncombo,nboot,njack);
      
      //load iboot
      int iboot_jack[100];
      FILE *fiboot=fopen(combine("%s/iboot",base_corrs_path[iens]).c_str(),"r");
      if(fiboot==NULL)
        {
          perror(combine("Error opening file iboot for ensamble %s",base_corrs_path[iens]).c_str());
          exit(1);
        }
      int nr=fread(iboot_jack,sizeof(int),100,fiboot);
      if(nr!=100)
        {
          perror(combine("Error loading iboot data for ensamble %s",base_corrs_path[iens]).c_str());
          exit(1);
        }
      
      //reading of data
      jvec tempm(ncombo,njack),tempz(ncombo,njack);
      tempm.load(M_path,0);
      
      //bootjacking
      for(int icombo=0;icombo<ncombo;icombo++)
	boot_from_jack(M[iens].data[icombo],tempm.data[icombo],iboot_jack);
    }
}
