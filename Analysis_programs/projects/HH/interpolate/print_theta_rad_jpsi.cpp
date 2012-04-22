#include "../HH_common.cpp"

int main()
{
  //load ensemble list path and data path
  FILE *an_input_file=open_file("analysis_pars","r");
  char ens_list_path[1024],base_M_path[1024];
  read_formatted_from_file_expecting(ens_list_path,an_input_file,"%s","ens_list_path");
  read_formatted_from_file_expecting(base_M_path,an_input_file,"%s","base_M_path");
  fclose(an_input_file);
  
  //load ensembles list and parameters
  char **base_corrs_path,**ens_name;
  int nens,*T,*ibeta,*nmass;
  double **mass,*sea_mass;
  load_ensembles_list(base_corrs_path,ens_name,nens,T,ibeta,nmass,mass,sea_mass,ens_list_path);
  
  //load all ensembles data
  bvec *aM1;
  bvec *aM2;
  
  load_all_ensembles_M(aM1,nens,T,ibeta,nmass,base_M_path,"VKVK",ens_name,base_corrs_path,"M");
  load_all_ensembles_M(aM2,nens,T,ibeta,nmass,base_M_path,"P5P5",ens_name,base_corrs_path,"M");

  for(int iens=0;iens<nens;iens++)
    {
      int ib=ibeta[iens];
      
      //prepare input
      boot q(nboot,njack);
      int im=1;
      {
	int ic=icombo(im,im,nmass[iens]);
	q=(aM1[iens][ic]*aM1[iens][ic]-aM2[iens][ic]*aM2[iens][ic])/(2*aM1[iens][ic]);
      }
      boot theta=T[iens]/2/M_PI/sqrt(3)*q;
      cout<<"iens: "<<iens<<", theta= "<<theta<<endl;
    }
  
  return 0;
}
