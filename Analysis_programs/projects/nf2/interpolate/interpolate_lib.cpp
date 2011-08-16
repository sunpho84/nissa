#include "../common.cpp"

//interpolate in the charm
bvec interpolate_charm(bvec vec,int nmass,int nlights,double *mass,int ibeta,const char *suff=NULL)
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

      if(suff!=NULL)
	{
	  ofstream outplot(combine("interpolating_charm_%s_l%d.xmg",suff,iml).c_str());
	  outplot<<"@type xydy"<<endl;
	  outplot<<"@s0 symbol 1"<<endl;
	  for(int imass=0;imass<nmass-nlights;imass++)
	    outplot<<temp_mass[imass]<<" "<<temp_vec[imass]<<endl;
	  outplot<<"&\n@type xydxdy"<<endl;
	  outplot<<amc_phys[ibeta].med()<<" "<<interpolated[iml].med()<<" ";
	  outplot<<amc_phys[ibeta].err()<<" "<<interpolated[iml].err()<<" ";
	  outplot.close();
	}
    }

  return interpolated;
}

//interpolate to the charm a degenerate meson
boot interpolate_charm_charm(bvec vec,int nmass,int nlights,double *mass,int ibeta,const char *suff=NULL)
{
  
  //temporary storage
  bvec temp_vec(nmass-nlights,nboot,njack);
  double temp_mass[nmass-nlights];
  //interpolate
  for(int imq=nlights;imq<nmass;imq++)
    {
      temp_vec[imq-nlights]=vec[icombo(imq,imq,nmass)];
      temp_mass[imq-nlights]=mass[imq];
    }
  
  return interpolate_single(temp_vec,temp_mass,amc_phys[ibeta]);
}

//interpolate in the strange
bvec interpolate_strange(bvec vec,int nmass,int nlights,double *mass,int ibeta)
{
  //output
  bvec interpolated(nlights,nboot,njack);
  
  //interpolate each light
  for(int iml=0;iml<nlights;iml++)
    {
      int nscart=min(nlights-2,iml+1);
      //temporary storage
      bvec temp_vec(nlights-nscart,nboot,njack);
      double temp_mass[nlights-nscart];
      //copy data for each light
      for(int imq=0;imq<nlights-nscart;imq++)
	{
	  temp_vec[imq]=vec[icombo(iml,imq+nscart,nmass)];
	  temp_mass[imq]=mass[imq+nscart];
	}
      interpolated[iml]=interpolate_single(temp_vec,temp_mass,ams_phys[ibeta]);
    }
  
  return interpolated;
}

//interpolate in the strange
void interpolate_many_strange(bvec *out,bvec in,int nmass,int nlights,double *mass,int ibeta,double *mint,int nint)
{
  for(int iref=0;iref<3;iref++)
    {
      //output
      out[iref]=bvec(nlights,nboot,njack);
      
      //interpolate each light
      for(int iml=0;iml<nlights;iml++)
	{
	  int nscart=min(nlights-2,iml+1);
	  //temporary storage
	  bvec temp_vec(nlights-nscart,nboot,njack);
	  double temp_mass[nlights-nscart];
	  //copy data for each light
	  for(int imq=0;imq<nlights-nscart;imq++)
	    {
	      temp_vec[imq]=in[icombo(iml,imq+nscart,nmass)];
	      temp_mass[imq]=mass[imq+nscart];
	    }
	  //cout<<iref<<" "<<mint[iref]<<endl;
	  out[iref][iml]=interpolate_single(temp_vec,temp_mass,mint[iref]*lat[ibeta]*Zp[ibeta]);
	}
    }
}

//interpolate in charm and in strange  
boot interpolate_charm_strange(bvec vec,int nmass,int nlights,double *mass,int ibeta,const char *suff=NULL)
{
  //first of all interpolate in the charm
  bvec charm_interpolated=interpolate_charm(vec,nmass,nlights,mass,ibeta,suff);
  
  //then interpolate in the strange
  boot out=interpolate_single(charm_interpolated,mass,ams_phys[ibeta]);
  
  if(suff!=NULL)
    {
      ofstream outplot(combine("interpolating_strange_%s.xmg",suff).c_str());
      outplot<<"@type xydy"<<endl;
      outplot<<"@s0 symbol 1"<<endl;
      for(int imass=0;imass<nlights;imass++)
	outplot<<mass[imass]<<" "<<charm_interpolated[imass]<<endl;
      outplot<<"&\n@type xydxdy"<<endl;
      outplot<<ams_phys[ibeta].med()<<" "<<out.med()<<" ";
      outplot<<ams_phys[ibeta].err()<<" "<<out.err()<<" ";
      outplot.close();
    }
  
  return out;
}

//load all data
void load_all_ensembles_MZ(bvec *&M,bvec *&Z,int &nens,int *&T,int *&ibeta,int *&nmass,const char *base_MZ_path,const char *obs_name,char **ens_name,char **base_corrs_path)
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
      
      //load iboot
      int iboot_jack[100];
      FILE *fiboot=fopen(combine("%s/iboot",base_corrs_path[iens]).c_str(),"r");
      if(fiboot==NULL)
	{
	  perror(combine("Error opening file iboot for ensamble %s",base_corrs_path[iens]).c_str());
	  exit(1);
	}
      int nr=fread(iboot_jack,sizeof(int),100,fiboot);
      if(nr!=100)
	{
	  perror(combine("Error loading iboot data for ensamble %s",base_corrs_path[iens]).c_str());
	  exit(1);
	}
      
      //reading of data
      jvec tempm(ncombo,njack),tempz(ncombo,njack);
      tempm.load(MZ_path,0);
      tempz.load(MZ_path,1);
      
      //bootjacking
      for(int icombo=0;icombo<ncombo;icombo++)
	{
	  boot_from_jack(M[iens].data[icombo],tempm.data[icombo],iboot_jack);
	  boot_from_jack(Z[iens].data[icombo],tempz.data[icombo],iboot_jack);
	}
    }
}
