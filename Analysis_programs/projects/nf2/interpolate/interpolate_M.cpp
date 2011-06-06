#include "../common.cpp"
#include "interpolate_lib.cpp"

int main()
{
  //load ensemble list path and data path
  FILE *an_input_file=open_file("analysis_pars","r");
  char ens_list_path[1024],base_MZ_path[1024],obs_name[1024],meson_name[1024];
  read_formatted_from_file_expecting(ens_list_path,an_input_file,"%s","ens_list_path");
  read_formatted_from_file_expecting(base_MZ_path,an_input_file,"%s","base_MZ_path");
  read_formatted_from_file_expecting(obs_name,an_input_file,"%s","obs_name");
  read_formatted_from_file_expecting(meson_name,an_input_file,"%s","meson_name");
  fclose(an_input_file);
  
  //load ensembles list and parameters
  char **base_corrs_path,**ens_name;
  int nens,*T,*ibeta,*iml_un,*nlights,*nmass;
  double **mass;
  load_ensembles_list(base_corrs_path,ens_name,nens,T,ibeta,nmass,mass,iml_un,nlights,ens_list_path);
  
  //load all ensembles data
  bvec *aM,*Z;
  load_all_ensembles_MZ(aM,Z,nens,T,ibeta,nmass,base_MZ_path,obs_name,ens_name,base_corrs_path);
  
  //compute M
  bvec M[nens];
  for(int iens=0;iens<nens;iens++)
    {
      int b=ibeta[iens];
      M[iens]=aM[iens]/lat[b];
    }
  
  bvec fint(nens,nboot,njack);
  for(int iens=0;iens<nens;iens++)
    {
      if(string(meson_name)==string("Ds"))
	fint.data[iens]=interpolate_charm_strange(M[iens],nmass[iens],nlights[iens],mass[iens],ibeta[iens]);
      
      if(string(meson_name)==string("D"))
	fint.data[iens]=interpolate_charm(M[iens],nmass[iens],nlights[iens],mass[iens],ibeta[iens])[iml_un[iens]];
      
      if(string(meson_name)==string("K"))
	fint.data[iens]=interpolate_strange(M[iens],nmass[iens],nlights[iens],mass[iens],ibeta[iens])[iml_un[iens]];
      
      if(string(meson_name)==string("Pi"))
	fint.data[iens]=M[iens][icombo(iml_un[iens],iml_un[iens],nmass[iens])];
      
      cout<<mass[iens][iml_un[iens]]<<" "<<fint.data[iens]<<endl;
    }
  
  if(string(meson_name)==string("Pi")) fint.write_to_binfile(combine("M_%s",meson_name).c_str());
  else fint.write_to_binfile(combine("interpolated_M_%s",meson_name).c_str());
  
  return 0;
}
