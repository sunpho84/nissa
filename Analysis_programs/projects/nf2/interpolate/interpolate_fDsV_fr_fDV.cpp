#include "../common.cpp"
#include "interpolate_lib.cpp"

int main()
{
  //load ensemble list path and data path
  FILE *an_input_file=open_file("analysis_pars","r");
  char ens_list_path[1024],base_MZ_path[1024],obs_name[1024];
  int mode;
  read_formatted_from_file_expecting(ens_list_path,an_input_file,"%s","ens_list_path");
  read_formatted_from_file_expecting(base_MZ_path,an_input_file,"%s","base_MZ_path");
  read_formatted_from_file_expecting((char*)&mode,an_input_file,"%d","mode");
  read_formatted_from_file_expecting(obs_name,an_input_file,"%s","obs_name");
  fclose(an_input_file);
  
  //load ensembles list and parameters
  char **base_corrs_path,**ens_name;
  int nens,*T,*ibeta,*iml_un,*nlights,*nmass;
  double **mass;
  load_ensembles_list(base_corrs_path,ens_name,nens,T,ibeta,nmass,mass,iml_un,nlights,ens_list_path);
  
  //load all ensembles data
  bvec *aM,*Z;
  load_all_ensembles_MZ(aM,Z,nens,T,ibeta,nlights,nmass,base_MZ_path,obs_name,ens_name,base_corrs_path,mode);
  
  cout<<"Interpolating for phi!"<<endl;
  
  //compute phi
  bvec phi[nens];
  for(int iens=0;iens<nens;iens++)
    {
      int b=ibeta[iens];
      
      phi[iens]=Za[b]*sqrt(Z[iens])/aM[iens]/lat[b]*sqrt(aM[iens]/lat[b]);
    }

  for(int iens=0;iens<nens;iens++)
    for(int ims=nmass[iens]-1;ims>=0;ims--)
      for(int imc=ims;imc<nmass[iens];imc++)
	{
	  int ic=icombo(imc,ims,nmass[iens],nlights[iens],mode);
	  phi[iens].data[ic]/=phi[iens].data[icombo(imc,iml_un[iens],nmass[iens],nlights[iens],mode)];
	}

  bvec phint(nens,nboot,njack);
  for(int iens=0;iens<nens;iens++)
    {
      phint.data[iens]=interpolate_charm_strange(phi[iens],nmass[iens],nlights[iens],mass[iens],ibeta[iens],mode);
      cout<<mass[iens][iml_un[iens]]<<" "<<phint.data[iens]<<endl;
    }
  
  phint.write_to_binfile("interpolated_phi_ratio_DsVDV");
  
  return 0;
}
