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
  bvec *aM_VECT,*Z_VECT;
  bvec *aM_PSEU,*Z_PSEU;
  load_all_ensembles_MZ(aM_VECT,Z_VECT,nens,T,ibeta,nmass,base_MZ_path,"VECT",ens_name);
  load_all_ensembles_MZ(aM_PSEU,Z_PSEU,nens,T,ibeta,nmass,base_MZ_path,"P5P5",ens_name);
  
  bvec Mint(nens,nboot,njack);
  for(int iens=0;iens<nens;iens++)
    {
      int ncombo=nmass[iens]*(nmass[iens]+1)/2;
      bvec X(ncombo,nboot,njack);
      for(int ims=0;ims<nmass[iens];ims++)
	{
	  //ofstream out(combine("%02d_%0.4f.xmg",iens,mass[iens][ims]).c_str());
	  
	  //prepare physical data
	  for(int imc=ims;imc<nmass[iens];imc++)
	    {
	      int ic=icombo(imc,ims,nmass[iens]);
	      X[ic]=(aM_VECT[iens][ic]-aM_PSEU[iens][ic])/aM_PSEU[iens][ic];
	      //out<<(mass[iens][imc]/lat[ibeta[iens]]/Zp[ibeta[iens]]).med()<<" "<<X[ic]<<endl;
	    }
	}
      
      Mint.data[iens]=interpolate_charm_strange(X,nmass[iens],nlights[iens],mass[iens],ibeta[iens]);
      
      cout<<(mass[iens][iml_un[iens]]/lat[ibeta[iens]]/Zp[ibeta[iens]]).med()<<" "<<Mint[iens]<<endl;
    }
  
  Mint.write_to_binfile("interpolated_M_Ds_ratdiff_V_PS");
  
  return 0;
}
