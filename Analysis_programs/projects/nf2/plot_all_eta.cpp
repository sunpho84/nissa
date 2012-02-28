#include "common.cpp"
#include "interpolate/interpolate_lib.cpp"

int main()
{
  //load ensemble list path and data path
  FILE *an_input_file=open_file("analysis_pars","r");
  int mode,iens;
  char ens_list_path[1024],base_MZ_path[1024],obs_name[1024];
  read_formatted_from_file_expecting(ens_list_path,an_input_file,"%s","ens_list_path");
  read_formatted_from_file_expecting(base_MZ_path,an_input_file,"%s","base_MZ_path");
  read_formatted_from_file_expecting((char*)&mode,an_input_file,"%d","mode");
  read_formatted_from_file_expecting((char*)&obs_name,an_input_file,"%s","obs_name");
  read_formatted_from_file_expecting((char*)&iens,an_input_file,"%d","iens");
  fclose(an_input_file);
  
  //load ensembles list and parameters
  char **base_corrs_path,**ens_name;
  int nens,*T,*ibeta,*iml_un,*nlights,*nmass;
  double **mass;
  load_ensembles_list(base_corrs_path,ens_name,nens,T,ibeta,nmass,mass,iml_un,nlights,ens_list_path);
  
  //load all ensembles data
  bvec *aM,*Z;
  load_all_ensembles_MZ(aM,Z,nens,T,ibeta,nlights,nmass,base_MZ_path,obs_name,ens_name,base_corrs_path,mode);
  
  ofstream out("Eta.xmg");
  out<<"@type xydy"<<endl;
  for(int im=0;im<nmass[iens];im++)
    out<<mass[iens][im]<<" "<<aM[iens][icombo(im,im,nmass[iens],nlights[iens],mode)]<<endl;
  out.close();
  
  return 0;
}
