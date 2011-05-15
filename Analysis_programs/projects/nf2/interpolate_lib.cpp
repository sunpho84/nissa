#include "common.cpp"

char data_list_file[1024];

void read_pars(const char *input)
{
  FILE *fin=open_file(input,"r");
  read_formatted_from_file_expecting(data_list_file,fin,"%s","data_list_file");
  
  fclose(fin);
}

int icombo(int im1,int im2)
{
  int ims=min(im1,im2);
  int imc=max(im1,im2);
  
  return ims*nmass-ims*(ims-1)/2+(imc-ims);
}

//interpolate in the strange or in the charm
bvec interpolate_strange_or_charm(bvec vec,boot amq)
{
  //output
  bvec interpolated(nlights,nboot,njack);
  
  //temporary storage
  bvec temp_vec(nmass,nboot,njack);
  
  //interpolate each light
  for(int iml=0;iml<nlights;iml++)
    {
      //copy data for each light
      for(int imq=0;imq<nmass;imq++) temp_vec[imq]=vec[icombo(iml,imq)];
      
      interpolated[iml]=interpolate_single(temp_vec,mass,amq);
    }
  
  return interpolated;
}

//interpolate in the charm
bvec interpolate_charm(bvec vec){return interpolate_charm_or_strange(vec,amc_phys[ibeta]);}

//interpolate in the strange
bvec interpolate_strange(bvec vec){return interpolate_charm_or_strange(vec,ams_phys[ibeta]);}
  
boot interpolate_charm_strange(bvec vec)
{
  //first of all interpolate in the charm
  bvec charm_interpolated=interpolate_charm(vec);
  
  //then interpolate in the strange
  return interpolate_single(charm_interpolated,mass,ams_phys[ibeta]);
}

int main()
{
  read_pars("input");
  read_input(data_list_file);
  init_latpars();
  
  bvec amD12(ncombo,nboot,njack);
  bvec afD12(ncombo,nboot,njack);
  
  //load all the corrs
  int ncombo=nmass*(nmass+1)/2;
  double *buf=new double[ncombo*T*(njack+1)];
  FILE *fin=open_file(combine("%s/%s",base_path,"P5P5").c_str(),"r");
  int stat=fread(buf,sizeof(double),ncombo_red*(njack+1)*T,fin);
  if(stat!=ncombo_red*(njack+1)*T)
    {
      cerr<<"Error loading data!"<<endl;
      exit(1);
    }
  
  int icombo=0;
  for(int ims=0;ims<nmass;ims++)
    for(int imc=ims;imc<nmass;imc++)
      {
	jvec corr(T,njack);
	corr.put(buf+icombo*T*(njack+1));
	
	int tmin=tmin_l;
	int tmax=tmax_l;
	if((ims>=1+nstrange)||(imc>=1+nstrange))
	  {
	    tmin=tmin_h;
	    tmax=tmax_h;
	  }
	
	jack tempm(njack),tempf(njack);
	fit(tempm,tempf,corr.simmetrized(1),combine("Ds_%d_%d",ims,imc).c_str(),tmin,tmax);

	boot_from_jack(amD12.data[icombo],tempm,ibeta);
	boot_from_jack(afD12.data[icombo],(mass[imc]+mass[ims])*tempf/(tempm*sinh(tempm)),ibeta);

	cout<<ims<<" "<<imc<<" "<<mass[ims]<<" "<<mass[imc]<<" "<<icombo<<" "<<amD12[icombo]<<" "<<tempf<<endl;
	
	icombo++;
      }
  
  boot mDs=interpolate_double(amD12)/lat[ibeta];
  boot fDs=interpolate_double(afD12)/lat[ibeta];
  
  mDs.write_to_binfile("results");
  fDs.append_to_binfile("results_fDs");

  cout<<mDs<<endl;
  cout<<fDs<<endl;

  return 0;
}
