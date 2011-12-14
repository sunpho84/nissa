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
  
  //prepare the list of mass
  int nh=nmass[2]-nlights[2];
  bvec mh(nh,nboot,njack);
  for(int ih=0;ih<nh;ih++)
    {
      mh[ih]=mass[1][ih+nlights[1]]/lat[1]/Zp[1];
      cout<<ih<<" "<<mh[ih]<<endl;
    }
  
  //compute M
  bvec M[nens];
  for(int iens=0;iens<nens;iens++)
    {
      int b=ibeta[iens];
      M[iens]=aM[iens]/lat[b];
    }
  
  
  bvec intM[nh];
  //prepare Mhl
  for(int ih=0;ih<nh;ih++)
    {
      intM[ih]=bvec(nens,nboot,njack);
      for(int iens=0;iens<nens;iens++)
	{
	  int ib=ibeta[iens];
	  int nl=nlights[ib];
	  int nm=nmass[ib];
	  int ni=nm-nl;
	  
	  //prepare input
	  double mi[ni];
	  bvec Mi(ni,nboot,njack);
	  for(int im=0;im<ni;im++)
	    {
	      mi[im]=mass[iens][im+nl];
	      Mi[im]=M[iens][im+nl];
	    } 
	  
	  //interpolate
	  boot x=mh[ih]*lat[ib]*Zp[ib];
	  intM[ih][iens]=interpolate_single(Mi,mi,x);

	  cout<<ih<<" iens="<<iens<<" x="<<x<<" y="<<intM[ih][iens]<<endl;
	}
      if(ih==0) intM[ih].write_to_binfile("intM");
      else      intM[ih].append_to_binfile("intM");
    }
  /*
//interpolate in the charm
bvec interpolate_many_charm(bvec in,int nmass,double *mass,int nint,double *mint,int ibeta)
{
  bvec out(nint,nboot,njack);
  for(int iint=0;iint<nint;iint++)
    out[iint]=
}

  */    
  
  return 0;
}
