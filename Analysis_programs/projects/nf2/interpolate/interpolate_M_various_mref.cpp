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
  
  double mint[3]={90.0/1000,95.0/1000,100.0/1000};
  bvec out(3*nens,nboot,njack);

  for(int iens=0;iens<nens;iens++)
    {
      bvec fint[3];
      if(string(meson_name)==string("K"))
	interpolate_many_strange(fint,sqr(M[iens]),nmass[iens],nlights[iens],mass[iens],ibeta[iens],mint,3);
      
      for(int iref=0;iref<3;iref++)
	{
	  //cout<<mass[iens][iml_un[iens]]<<" "<<fint[iref][iml_un[iens]]<<endl;
	  out[iref*nens+iens]=fint[iref][iml_un[iens]];
	}
    }
  out.write_to_binfile(combine("interpolated_M2_%s",meson_name).c_str());
  
  return 0;
}
