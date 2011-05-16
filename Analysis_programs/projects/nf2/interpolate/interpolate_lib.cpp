#include "../common.cpp"

//interpolate in the charm
bvec interpolate_charm(bvec vec,int nmass,int nlights,double *mass,int ibeta)
{
  //output
  bvec interpolated(nlights,nboot,njack);
  
  //temporary storage
  bvec temp_vec(nmass-nlights,nboot,njack);
  double temp_mass[nmass-nlights];
  //interpolate each light
  for(int iml=0;iml<nlights;iml++)
    {
      //copy data for each light
      for(int imq=nlights;imq<nmass;imq++)
	{
	  temp_vec[imq-nlights]=vec[icombo(iml,imq,nmass)];
	  temp_mass[imq-nlights]=mass[imq];
	}
      
      interpolated[iml]=interpolate_single(temp_vec,temp_mass,amc_phys[ibeta]);
    }
  
  return interpolated;
}

//interpolate in the strange
bvec interpolate_strange(bvec vec,int nmass,int nlights,double *mass,int ibeta)
{
  //output
  bvec interpolated(nlights,nboot,njack);
  
  //temporary storage
  bvec temp_vec(nmass-nlights,nboot,njack);
  double temp_mass[nmass-nlights];
  //interpolate each light
  for(int iml=0;iml<nlights;iml++)
    {
      //copy data for each light
      for(int imq=0;imq<nlights;imq++)
	{
	  temp_vec[imq]=vec[icombo(iml,imq,nmass)];
	  temp_mass[imq]=mass[imq];
	}
      
      interpolated[iml]=interpolate_single(temp_vec,temp_mass,ams_phys[ibeta]);
    }
  
  return interpolated;
}

//interpolate in charm and in strange  
boot interpolate_charm_strange(bvec vec,int nmass,int nlights,double *mass,int ibeta)
{
  //first of all interpolate in the charm
  bvec charm_interpolated=interpolate_charm(vec,nmass,nlights,mass,ibeta);
  
  //then interpolate in the strange
  return interpolate_single(charm_interpolated,mass,ams_phys[ibeta]);
}

//load all data
void load_all_ensambles_MZ(const char *ens_list_path,const char *obs_name,bvec *&M,bvec *&Z,int &nens,int *&T,int *&ibeta,int *&nmass,double **&mass,int *&iml_un,int *&nlights)
{
  init_latpars();
  
  FILE *ens_list_file=open_file(ens_list_path,"r");
  read_formatted_from_file_expecting((char*)(&nens),ens_list_file,"%d","nens");
  
  char ens_path_base[1024],data_path_base[1024];
  read_formatted_from_file_expecting((char*)(&ens_path_base),ens_list_file,"%s","ens_path_base");
  read_formatted_from_file_expecting((char*)(&data_path_base),ens_list_file,"%s","data_path_base");
  
  M=(bvec*)malloc(nens*sizeof(bvec));
  Z=(bvec*)malloc(nens*sizeof(bvec));
  T=(int*)malloc(nens*sizeof(int));
  ibeta=(int*)malloc(nens*sizeof(int));
  nmass=(int*)malloc(nens*sizeof(int));
  mass=(double**)malloc(nens*sizeof(double*));
  iml_un=(int*)malloc(nens*sizeof(int));
  nlights=(int*)malloc(nens*sizeof(int));
  
  for(int iens=0;iens<nens;iens++)
    {
      char ens_path_pars[1024],base_path[1024];
      read_formatted_from_file(ens_path_pars,ens_list_file,"%s","ens_path_pars");
      
      char ens_path[1024],data_path[1024];
      sprintf(ens_path,"%s/%s",ens_path_base,ens_path_pars);
      sprintf(data_path,"%s/%s/%s/results",data_path_base,obs_name,ens_path_pars);

      cout<<"Reading ensamble: "<<ens_path<<endl;
      read_ensamble_pars(base_path,T[iens],ibeta[iens],nmass[iens],mass[iens],iml_un[iens],nlights[iens],combine("%s/data_list",ens_path).c_str());

      int ncombo=nmass[iens]*(nmass[iens]+1)/2;
      M[iens].create(ncombo,nboot,njack);
      Z[iens].create(ncombo,nboot,njack);
      
      //reading of data
      jvec tempm(ncombo,njack);
      jvec tempz(ncombo,njack);
      tempm.load(data_path,0);
      tempz.load(data_path,1);
      
      //bootjacking
      for(int icombo=0;icombo<ncombo;icombo++)
	{
	  boot_from_jack(M[iens].data[icombo],tempm.data[icombo],ibeta[iens]);
	  boot_from_jack(Z[iens].data[icombo],tempz.data[icombo],ibeta[iens]);
	}
    }
  
  fclose(ens_list_file);
}
