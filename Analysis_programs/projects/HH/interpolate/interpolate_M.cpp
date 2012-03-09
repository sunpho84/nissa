#include "../HH_common.cpp"

int main()
{
  //load ensemble list path and data path
  FILE *an_input_file=open_file("analysis_pars","r");
  char ens_list_path[1024],base_M_path[1024],obs_name[1024];
  read_formatted_from_file_expecting(ens_list_path,an_input_file,"%s","ens_list_path");
  read_formatted_from_file_expecting(base_M_path,an_input_file,"%s","base_M_path");
  read_formatted_from_file_expecting(obs_name,an_input_file,"%s","obs_name");
  fclose(an_input_file);
  
  //load ensembles list and parameters
  char **base_corrs_path,**ens_name;
  int nens,*T,*ibeta,*nmass;
  double **mass,*sea_mass;
  load_ensembles_list(base_corrs_path,ens_name,nens,T,ibeta,nmass,mass,sea_mass,ens_list_path);
  
  //load all ensembles data
  bvec *aM;
  if(strlen(obs_name)==4)
    load_all_ensembles_M(aM,nens,T,ibeta,nmass,base_M_path,obs_name,ens_name,base_corrs_path,"M");
  else
    {
      char sep[2]="_";
      char *last;
      char *first=strtok_r(obs_name,sep,&last);
      char *opera=strtok_r(NULL,sep,&last);
      char *secon=strtok_r(NULL,sep,&last);
      cout<<first<<" "<<opera<<" "<<secon<<endl;
      bvec *aM2;
      load_all_ensembles_M(aM,nens,T,ibeta,nmass,base_M_path,first,ens_name,base_corrs_path,"M");
      load_all_ensembles_M(aM2,nens,T,ibeta,nmass,base_M_path,secon,ens_name,base_corrs_path,"M");
      for(int iens=0;iens<nens;iens++)
	if(string(opera)=="minus")
	  aM[iens]-=aM2[iens];
    }
  
  init_latpars();
  
  bvec Mout(nref_hmass*nens,nboot,njack);
  //prepare Mhl
  for(int iens=0;iens<nens;iens++)
    {
      int ib=ibeta[iens];
      
      //prepare input
      bvec mass_out(nref_hmass,nboot,njack);
      bvec Min(nmass[iens],nboot,njack);
      for(int im=0;im<nmass[iens];im++)
	//if(string("VKVK_minus_P5P5")!=obs_name)
	  Min[im]=aM[iens][icombo(im,im,nmass[iens])]/lat[ib];
      //else
      //Min[im]=aM[iens][icombo(im,im,nmass[iens])];
      for(int iref_hmass=0;iref_hmass<nref_hmass;iref_hmass++)
	mass_out[iref_hmass]=ref_hmass[iref_hmass]*lat[ib]*Zp[ib];
      
      //interpolate
      bvec temp_Mout=interpolate_multi(mass[iens],Min,mass_out,combine("interpolating_ens%02d.xmg",iens).c_str());
      for(int iref_hmass=0;iref_hmass<nref_hmass;iref_hmass++)
	{
	  Mout[iens+iref_hmass*nens]=temp_Mout[iref_hmass];
	  cout<<iref_hmass<<" iens="<<iens<<" x="<<mass_out[iref_hmass]<<" y="<<temp_Mout[iref_hmass]<<endl;
	}
    }
  
  ofstream out390("390charm.xmg");
  out390<<"@type xydy"<<endl;
  for(int iens=0;iens<nens;iens++)
    if(ibeta[iens]==1)
      out390<<sea_mass[iens]<<" "<<aM[iens][icombo(1,1,nmass[iens])]<<endl;
    
  Mout.write_to_binfile("intM");
  
  return 0;
}
