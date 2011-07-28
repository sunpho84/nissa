#include "../common.cpp"
#include "interpolate_lib.cpp"

int main()
{
  //load ensemble list path and data path
  FILE *an_input_file=open_file("analysis_pars","r");
  char ens_list_path[1024],base_MZ_path[1024];
  read_formatted_from_file_expecting(ens_list_path,an_input_file,"%s","ens_list_path");
  read_formatted_from_file_expecting(base_MZ_path,an_input_file,"%s","base_MZ_path");
  fclose(an_input_file);
  
  //load ensembles list and parameters
  char **base_corrs_path,**ens_name;
  int nens,*T,*ibeta,*iml_un,*nlights,*nmass;
  double **mass;
  load_ensembles_list(base_corrs_path,ens_name,nens,T,ibeta,nmass,mass,iml_un,nlights,ens_list_path);
  
  //load all ensembles data
  bvec *aM_PS,*Z_PS;
  load_all_ensembles_MZ(aM_PS,Z_PS,nens,T,ibeta,nmass,base_MZ_path,"P5P5",ens_name,base_corrs_path);
  
  bvec fint(nens,nboot,njack);
  
  double ph_EtaC=2980.3;
  double ph_Ds=1968.5;
  double ph_ratdiff=(ph_EtaC/2-ph_Ds);
  
  for(int iens=0;iens<nens;iens++)
    {
      boot Ds=interpolate_charm_strange(aM_PS[iens],nmass[iens],nlights[iens],mass[iens],ibeta[iens]);
      boot EtaC=interpolate_charm_charm(aM_PS[iens],nmass[iens],nlights[iens],mass[iens],ibeta[iens]);
      
      fint.data[iens]=(EtaC/2-Ds)/lat[ibeta[iens]];
      
      cout<<mass[iens][iml_un[iens]]<<" "<<fint.data[iens]<<endl;
    }
  
  fint.write_to_binfile("Mass_ratdiff.dat");
  
  cout<<ph_ratdiff<<endl;
  
  return 0;
}
