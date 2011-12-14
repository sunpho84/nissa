#include "../common.cpp"
#include "interpolate_lib.cpp"

//load all data
void load_all_ensembles_X(bvec *&X,int &nens,int *&T,int *&ibeta,int *&nlights,int *&nmass,const char *base_X_path,char **ens_name,char **base_corrs_path,int mode)
{
  init_latpars();
  
  //allocate room for m and z
  X=(bvec*)malloc(nens*sizeof(bvec));

  //Loop over ensembles. Data is supposed to be stored in a file named [base_X_path]/TFRV/[ens_name]/results
  for(int iens=0;iens<nens;iens++)
    {
      char X_path[1024];
      sprintf(X_path,"%s/TFRV/%s/results",base_X_path,ens_name[iens]);

      cout<<"Reading ensemble: "<<ens_name[iens]<<" from file: "<<X_path<<endl;

      //allocate room 
      int ncombo;
      switch(mode)
        {
        case 0:
          ncombo=nmass[iens]*(nmass[iens]+1)/2;
          break;
        case 1:
          ncombo=nlights[iens]*nmass[iens];
          break;
        case 2:
          ncombo=nlights[iens]*(nmass[iens]-(nlights[iens]-1)/2);
          break;
        default:
          cerr<<"Error, unkwnown mode!"<<endl;
          ncombo=0;
          break;
        }
      
      X[iens].create(ncombo,nboot,njack);
      
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
      tempm.load(X_path,0);
      
      //bootjacking
      for(int icombo=0;icombo<ncombo;icombo++)
	boot_from_jack(X[iens].data[icombo],tempm.data[icombo],iboot_jack);
    }
}

int main()
{
  //load ensemble list path and data path
  FILE *an_input_file=open_file("analysis_pars","r");
  int mode;
  char ens_list_path[1024],base_X_path[1024],meson_name[1024];
  read_formatted_from_file_expecting(ens_list_path,an_input_file,"%s","ens_list_path");
  read_formatted_from_file_expecting(base_X_path,an_input_file,"%s","base_X_path");
  read_formatted_from_file_expecting((char*)&mode,an_input_file,"%d","mode");
  read_formatted_from_file_expecting(meson_name,an_input_file,"%s","meson_name");
  fclose(an_input_file);
  
  //load ensembles list and parameters
  char **base_corrs_path,**ens_name;
  int nens,*T,*ibeta,*iml_un,*nlights,*nmass;
  double **mass;
  load_ensembles_list(base_corrs_path,ens_name,nens,T,ibeta,nmass,mass,iml_un,nlights,ens_list_path);
  
  //load all ensembles data
  bvec *X;
  load_all_ensembles_X(X,nens,T,ibeta,nlights,nmass,base_X_path,ens_name,base_corrs_path,mode);
  
  //compute f
  bvec x[nens];
  for(int iens=0;iens<nens;iens++)
    {
      int b=ibeta[iens];
      
      x[iens]=X[iens];
    }
  
  bvec fint(nens,nboot,njack);
  for(int iens=0;iens<nens;iens++)
    {
      if(string(meson_name)==string("Ds"))
	fint.data[iens]=interpolate_charm_strange(x[iens],nmass[iens],nlights[iens],mass[iens],ibeta[iens],mode,combine("out%02d",iens).c_str());
      
      if(string(meson_name)==string("D"))
	fint.data[iens]=interpolate_charm(x[iens],nmass[iens],nlights[iens],mass[iens],ibeta[iens],mode,combine("out%02d",iens).c_str())[iml_un[iens]];
      
      if(string(meson_name)==string("K"))
	fint.data[iens]=interpolate_unitary_light_strange(x[iens],nmass[iens],nlights[iens],iml_un[iens],mass[iens],ibeta[iens],mode,combine("out%02d",iens).c_str());
      
      if(string(meson_name)==string("Pi"))
	fint.data[iens]=x[iens][icombo(iml_un[iens],iml_un[iens],nmass[iens],nlights[iens],mode)];
      
      cout<<iens<<" "<<ibeta[iens]<<" "<<mass[iens][iml_un[iens]]<<" "<<fint.data[iens]<<endl;
    }
  
  fint.write_to_binfile(combine("interpolated_X_%s",meson_name).c_str());
  
  return 0;
}
