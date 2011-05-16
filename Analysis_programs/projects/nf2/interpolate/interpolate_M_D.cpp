#include "../common.cpp"
#include "interpolate_lib.cpp"

int main()
{
  bvec *aM,*Z;
  int nens,*T,*ibeta,*iml_un,*nlights,*nmass;
  double **mass;
  
  //load all ensembles parameters
  char **base_corrs_paths;
  load_ensemble_list(base_corrs_paths,ens_name,nens,T,ibeta,nmass,mass,iml_un,nlights,ens_list_path);
  

  for(int iens=0;iens<nens;iens++)
    {
      int ncombo=nmass[iens]*(nmass[iens]+1)/2;
      bvec M(ncombo,nboot,njack);
      for(int ims=0;ims<nmass[iens];ims++)
	{
	  ofstream out(combine("%02d_%0.4f.xmg",iens,mass[iens][ims]).c_str());
	  for(int imc=ims;imc<nmass[iens];imc++)
	    {
	      int ic=icombo(imc,ims,nmass[iens]);
	      M[ic]=aM[iens][ic]/lat[ibeta[iens]];
	      out<<(mass[iens][imc]/lat[ibeta[iens]]/Zp[ibeta[iens]]).med()<<" "<<M[ic]<<endl;
	    }
	  cout<<(mass[iens][ims]/lat[ibeta[iens]]/Zp[ibeta[iens]]).med()<<" "<<interpolate_charm(M,nmass[iens],nlights[iens],mass[iens],ibeta[iens])<<endl;
	}
    }
  
    return 0;
}
