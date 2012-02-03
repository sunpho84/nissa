#include "../common.cpp"
#include "interpolate_lib.cpp"

int main()
{
  //load ensemble list path and data path
  FILE *an_input_file=open_file("analysis_pars","r");
  int mode;
  char ens_list_path[1024],base_MZ_path[1024],obs_name[1024],meson_name[1024];
  read_formatted_from_file_expecting(ens_list_path,an_input_file,"%s","ens_list_path");
  read_formatted_from_file_expecting(base_MZ_path,an_input_file,"%s","base_MZ_path");
  read_formatted_from_file_expecting((char*)&mode,an_input_file,"%d","mode");
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
  load_all_ensembles_MZ(aM,Z,nens,T,ibeta,nlights,nmass,base_MZ_path,obs_name,ens_name,base_corrs_path,mode);
  
  //compute f
  bvec M[nens];
  bvec f[nens];
  bvec phi[nens];
  bvec x[nens];
  for(int iens=0;iens<nens;iens++)
    {
      int b=ibeta[iens];
      
      M[iens]=aM[iens]/lat[b];
      
      if(string(obs_name)==string("P5P5"))
	{
	  f[iens]=sqrt(Z[iens])/(sinh(aM[iens])*aM[iens])/lat[b];
	  phi[iens]=f[iens]*sqrt(M[iens]);
	  int lim1;
	  if(mode==0) lim1=nmass[iens];
	  else lim1=nlights[iens];
	  for(int ims=0;ims<lim1;ims++)
	    for(int imc=ims;imc<nmass[iens];imc++)
	      {
		int ic=icombo(ims,imc,nmass[iens],nlights[iens],mode);
		f[iens].data[ic]*=mass[iens][ims]+mass[iens][imc];
		phi[iens].data[ic]*=mass[iens][ims]+mass[iens][imc];
		cout<<iens<<" ens, ic="<<ic<<", il="<<mass[iens][ims]<<", ih="<<mass[iens][imc]<<", Z="<<Z[iens][ic]<<", M="<<M[iens][ic]<<endl;
	      }
	    x[iens]=phi[iens];
	}
      
      if(string(obs_name)==string("S0S0"))
	{
	  f[iens]=(1/Zp_fr_Zs[b])*sqrt(Z[iens])/(sinh(aM[iens])*aM[iens])/lat[b];
	  phi[iens]=f[iens]*sqrt(M[iens]);
	  int lim1;
	  if(mode==0) lim1=nmass[iens];
	  else lim1=nlights[iens];
	  for(int ims=0;ims<lim1;ims++)
	    for(int imc=ims;imc<nmass[iens];imc++)
	      {
		int ic=icombo(ims,imc,nmass[iens],nlights[iens],mode);
		f[iens].data[ic]*=mass[iens][imc]-mass[iens][ims];
		phi[iens].data[ic]*=mass[iens][imc]-mass[iens][ims];
		cout<<iens<<" ens, ic="<<ic<<", il="<<mass[iens][ims]<<", ih="<<mass[iens][imc]<<", Z="<<Z[iens][ic]<<", M="<<M[iens][ic]<<endl;
	      }
	  x[iens]=phi[iens];
	}
      
      if(string(obs_name)==string("VKVK"))
	{
	  x[iens]=f[iens]=Za[b]*sqrt(Z[iens])/aM[iens]/lat[b];
	  x[iens]*=sqrt(aM[iens]/lat[b]);
	}
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
  
  if(string(meson_name)==string("K")) fint.write_to_binfile(combine("interpolated_f_%s",meson_name).c_str());
  else
    if(string(meson_name)==string("Pi")) fint.write_to_binfile(combine("f_%s",meson_name).c_str());
    else fint.write_to_binfile(combine("interpolated_phi_%s",meson_name).c_str());
  
  return 0;
}
