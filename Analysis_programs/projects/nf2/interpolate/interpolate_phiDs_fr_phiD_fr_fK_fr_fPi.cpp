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
  load_all_ensembles_MZ(aM,Z,nens,T,ibeta,nmass,base_MZ_path,obs_name,ens_name,base_corrs_path);
  
  //compute f and phi
  bvec f[nens],phi[nens],ratio[nens];
  for(int iens=0;iens<nens;iens++)
    {
      int b=ibeta[iens];
      f[iens]=sqrt(Z[iens])/(sinh(aM[iens])*aM[iens])/lat[b];
      phi[iens]=f[iens]*sqrt(aM[iens]/lat[b]);
      for(int ims=0;ims<nmass[iens];ims++)
	for(int imc=ims;imc<nmass[iens];imc++)
	  {
	    int ic=icombo(imc,ims,nmass[iens]);
	    phi[iens][ic]*=mass[iens][ims]+mass[iens][imc];
	    f[iens][ic]*=mass[iens][ims]+mass[iens][imc];
	  }

      ratio[iens]=phi[iens];
      for(int ims=0;ims<nmass[iens];ims++)
	for(int imc=ims;imc<nmass[iens];imc++)
	  {
	    int nm=nmass[iens];
	    boot phiD=phi[iens][icombo(iml_un[iens],imc,nm)];
	    boot phiDs=phi[iens][icombo(ims,imc,nm)];
	    boot fPi=f[iens][icombo(iml_un[iens],iml_un[iens],nm)];
	    boot fK=f[iens][icombo(iml_un[iens],ims,nm)];
	    
	    ratio[iens][icombo(ims,imc,nm)]=(phiDs/phiD)/(fK/fPi);
	  }
    }

  bvec ratio_int(nens,nboot,njack);
  for(int iens=0;iens<nens;iens++)
    ratio_int[iens]=interpolate_charm_strange(ratio[iens],nmass[iens],nlights[iens],mass[iens],ibeta[iens],combine("ens%d",iens).c_str());

  ratio_int.write_to_binfile("interpolated_ratio");
  
  return 0;
}
