#include "../../nf2/common_pars.cpp"
#include "../../nf2/common.cpp"
#include "../../nf2/interpolate/interpolate_lib.cpp"

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
  
  init_latpars();
  
  //compute M
  bvec M[nens];
  for(int iens=0;iens<nens;iens++)
    {
      int b=ibeta[iens];
      M[iens]=aM[iens]/lat[b];
    }
  
  bvec intM(nref_hmass*nens,nboot,njack);
  //prepare Mhl
  for(int iens=0;iens<nens;iens++)
    {
      int ib=ibeta[iens];
      int nl=nlights[iens];
      int nm=nmass[iens];
      int nin=nm-nl;
      
      //prepare input
      double xin[nin];
      bvec yin(nin,nboot,njack);
      for(int im=0;im<nin;im++)
	{
	  xin[im]=mass[iens][im+nl];
	  yin[im]=M[iens][im+nl];
	} 
      
      //prepare output
      bvec xout(nref_hmass,nboot,njack),yout(nref_hmass,nboot,njack);
      for(int iref_hmass=0;iref_hmass<nref_hmass;iref_hmass++) xout[iref_hmass]=ref_hmass[iref_hmass]*lat[ib]*Zp[ib];
      
      //interpolate
      yout=interpolate_multi(xin,yin,xout,combine("interpolating_ens%02d.xmg",iens).c_str());
      for(int iref_hmass=0;iref_hmass<nref_hmass;iref_hmass++)
	{
	  intM[iens+iref_hmass*nens]=yout[iref_hmass];
	  cout<<iref_hmass<<" iens="<<iens<<" x="<<xout[iref_hmass].med()<<" y="<<yout[iref_hmass]<<endl;
	}
    }
  
  intM.write_to_binfile("intM");
  
  return 0;
}
