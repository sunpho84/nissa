#include "../common.cpp"
#include "interpolate_lib.cpp"

int main()
{
  //load ensemble list path and data path
  FILE *an_input_file=open_file("analysis_pars","r");
  char ens_list_path[1024],base_MZ_path[1024],obs_name[1024];
  read_formatted_from_file_expecting(ens_list_path,an_input_file,"%s","ens_list_path");
  read_formatted_from_file_expecting(base_MZ_path,an_input_file,"%s","base_MZ_path");
  read_formatted_from_file_expecting(obs_name,an_input_file,"%s","obs_name");
  fclose(an_input_file);
  
  //load ensembles list and parameters
  char **base_corrs_path,**ens_name;
  int nens,*T,*ibeta,*iml_un,*nlights,*nmass;
  double **mass;
  load_ensembles_list(base_corrs_path,ens_name,nens,T,ibeta,nmass,mass,iml_un,nlights,ens_list_path);
  
  //load all ensembles data
  bvec *aM,*Z;
  load_all_ensembles_MZ(aM,Z,nens,T,ibeta,nmass,base_MZ_path,obs_name,ens_name);
  
  //compute phi
  bvec phiD_s[nens];
  for(int iens=0;iens<nens;iens++)
    {
      int b=ibeta[iens];
      phiD_s[iens]=sqrt(Z[iens])/(sinh(aM[iens])*aM[iens])/lat[b]*sqrt(aM[iens]/lat[b]);
      for(int ims=0;ims<nmass[iens];ims++)
	for(int imc=ims;imc<nmass[iens];imc++)
	  {
	    int ic=icombo(imc,ims,nmass[iens]);
	    phiD_s[iens].data[ic]*=mass[iens][ims]+mass[iens][imc];
	  }
    }
  
  bvec phiD_s_int(nens,nboot,njack);
  for(int iens=0;iens<nens;iens++)
    {
      phiD_s_int.data[iens]=interpolate_charm_strange(phiD_s[iens],nmass[iens],nlights[iens],mass[iens],ibeta[iens]);
      
      cout<<(mass[iens][iml_un[iens]]/lat[ibeta[iens]]/Zp[ibeta[iens]]).med()<<" "<<phiD_s_int[iens]<<endl;
    }
  
  phiD_s_int.write_to_binfile("interpolated_phiD_s");
  
  return 0;
}
