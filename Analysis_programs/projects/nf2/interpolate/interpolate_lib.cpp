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
      
      interpolated[iml]=interpolate_single(temp_vec,temp_mass,mc_phys*lat_med[ibeta]*Zp_med[ibeta]);
    }
  
  return interpolated;
}

//interpolate in the strange
bvec interpolate_strange(bvec vec,int nmass,int nlights,double *mass,int ibeta)
{
  //output
  bvec interpolated(nlights,nboot,njack);
  
  //temporary storage
  bvec temp_vec(nlights,nboot,njack);
  double temp_mass[nlights];
  //interpolate each light
  for(int iml=0;iml<nlights;iml++)
    {
      //copy data for each light
      for(int imq=0;imq<nlights;imq++)
	{
	  temp_vec[imq]=vec[icombo(iml,imq,nmass)];
	  temp_mass[imq]=mass[imq];
	}
      
      interpolated[iml]=interpolate_single(temp_vec,temp_mass,ms_phys*lat_med[ibeta]*Zp_med[ibeta]);
    }
  
  return interpolated;
}

//interpolate in charm and in strange  
boot interpolate_charm_strange(bvec vec,int nmass,int nlights,double *mass,int ibeta)
{
  //first of all interpolate in the charm
  bvec charm_interpolated=interpolate_charm(vec,nmass,nlights,mass,ibeta);
  
  //then interpolate in the strange
  return interpolate_single(charm_interpolated,mass,ms_phys*lat_med[ibeta]*Zp_med[ibeta]);
}

//load all data
void load_all_ensembles_MZ(bvec *&M,bvec *&Z,int &nens,int *&T,int *&ibeta,int *&nmass,const char *base_MZ_path,const char *obs_name,char **ens_name)
{
  init_latpars();
  
  //allocate room for m and z
  M=(bvec*)malloc(nens*sizeof(bvec));
  Z=(bvec*)malloc(nens*sizeof(bvec));
  
  //Loop over ensembles. Data is supposed to be stored in a file named [base_MZ_path]/[obs_name]/[ens_name]/results
  for(int iens=0;iens<nens;iens++)
    {
      char MZ_path[1024];
      sprintf(MZ_path,"%s/%s/%s/results",base_MZ_path,obs_name,ens_name[iens]);
      
      cout<<"Reading ensemble: "<<ens_name[iens]<<" from file: "<<MZ_path<<endl;
      
      //allocate room or M and Z2
      int ncombo=nmass[iens]*(nmass[iens]+1)/2;
      M[iens].create(ncombo,nboot,njack);
      Z[iens].create(ncombo,nboot,njack);
      
      //reading of data
      jvec tempm(ncombo,njack),tempz(ncombo,njack);
      tempm.load(MZ_path,0);
      tempz.load(MZ_path,1);
      
      //bootjacking
      for(int icombo=0;icombo<ncombo;icombo++)
	{
	  boot_from_jack(M[iens].data[icombo],tempm.data[icombo],ibeta[iens]);
	  boot_from_jack(Z[iens].data[icombo],tempz.data[icombo],ibeta[iens]);
	}
    }
}
