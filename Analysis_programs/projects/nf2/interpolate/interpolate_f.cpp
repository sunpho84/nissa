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
  load_all_ensembles_MZ(aM,Z,nens,T,ibeta,nmass,base_MZ_path,obs_name,ens_name);
  
  //compute f
  bvec f[nens];
  for(int iens=0;iens<nens;iens++)
    {
      int b=ibeta[iens];
      
      if(string(obs_name)==string("P5P5"))
	{
	  f[iens]=sqrt(Z[iens])/(sinh(aM[iens])*aM[iens])/lat[b];
	  for(int ims=0;ims<nmass[iens];ims++)
	    for(int imc=ims;imc<nmass[iens];imc++)
	      {
		int ic=icombo(imc,ims,nmass[iens]);
		f[iens].data[ic]*=mass[iens][ims]+mass[iens][imc];
	      }
	}
      double r0=0.450/hc;
      
      if(string(obs_name)==string("VECT")) f[iens]=Za[b]*sqrt(Z[iens])/aM[iens]/lat[b];
    }
  
  bvec fint(nens,nboot,njack);
  for(int iens=0;iens<nens;iens++)
    {
      if(string(meson_name)==string("Ds"))
	fint.data[iens]=interpolate_charm_strange(f[iens],nmass[iens],nlights[iens],mass[iens],ibeta[iens]);
      
      if(string(meson_name)==string("D"))
	fint.data[iens]=interpolate_charm(f[iens],nmass[iens],nlights[iens],mass[iens],ibeta[iens])[iml_un[iens]];
	  
      if(string(meson_name)==string("K"))
	fint.data[iens]=interpolate_strange(f[iens],nmass[iens],nlights[iens],mass[iens],ibeta[iens])[iml_un[iens]];
      
      cout<<mass[iens][iml_un[iens]]<<" "<<fint.data[iens]<<endl;
    }
  
  fint.write_to_binfile(combine("interpolated_f_%s",meson_name).c_str());
  
  return 0;
}
