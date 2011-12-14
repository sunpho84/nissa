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

  //output
  bvec rat_shint[nh];
  for(int ih=0;ih<nh;ih++) rat_shint[ih]=bvec(nens,nboot,njack);
  //prepare the y for the interpolation
  bvec af[nens],aphi[nens];
  for(int iens=0;iens<nens;iens++)
    {
      //int b=ibeta[iens];
      af[iens]=sqrt(Z[iens])/(aM[iens]*sinh(aM[iens]));
      aphi[iens]=af[iens]*sqrt(aM[iens]);
      for(int iml=0;iml<nlights[iens];iml++)
	for(int imh=0;imh<nmass[iens];imh++)
	  {
	    int ic=iml*nmass[iens]+imh;
	    af[iens][ic]*=(mass[iens][iml]+mass[iens][imh]);
	    aphi[iens][ic]*=(mass[iens][iml]+mass[iens][imh]);
	    cout<<iens<<" "<<iml<<" "<<imh<<" f: "<<af[iens][ic]/lat[ibeta[iens]]<<" phi: "<<aphi[iens][ic]/pow(lat[ibeta[iens]],1.5)<<endl;
	  }
    }
  
  for(int iens=0;iens<nens;iens++)
    {
      //now interpolate to the strange each heavy
      double *xh=mass[iens]+nlights[iens];
      bvec rat_sint(nmass[iens]-nlights[iens],nboot,njack);

      ofstream outh(combine("hint_%d.xmg",iens).c_str());
      outh<<"@type xydy"<<endl;

      for(int ih=nlights[iens];ih<nmass[iens];ih++)
	{
	  double xto[nlights[iens]-1];
	  bvec yto(nlights[iens]-1,nboot,njack);
	  
	  ofstream out(combine("stint_%d_%d.xmg",iens,ih).c_str());
	  out<<"@type xydy"<<endl;
	  for(int iml=1;iml<nlights[iens];iml++)
	    {
	      xto[iml-1]=mass[iens][iml];
	      yto[iml-1]=aphi[iens][iml*nmass[iens]+ih]/pow(lat[ibeta[iens]],1.5);
	      out<<mass[iens][iml]<<" "<<yto[iml-1]<<endl;
	    }
	  out<<"&"<<endl;
	  rat_sint.data[ih-nlights[iens]]=interpolate_single(yto,xto,ams_phys[ibeta[iens]]);
	  out<<ams_phys[ibeta[iens]].med()<<" "<<rat_sint.data[ih-nlights[iens]]<<endl;

	  outh<<mass[iens][ih]<<" "<<rat_sint[ih-nlights[iens]]<<endl;
	}
      
      outh<<"&"<<endl;
      
      for(int ih=0;ih<nh;ih++)
	{
	  boot x=mh[ih]*lat[ibeta[iens]]*Zp[ibeta[iens]];
	  rat_shint[ih][iens]=interpolate_single(rat_sint,xh,x);
	  
	  outh<<x.med()<<" "<<rat_shint[ih][iens]<<endl;
	}
    }
  
  for(int ih=0;ih<nh;ih++)
    if(ih==0) rat_shint[ih].write_to_binfile("results");
    else      rat_shint[ih].append_to_binfile("results");

  return 0;
}
