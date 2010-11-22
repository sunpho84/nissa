#include <vector>
#include <utility>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cassert>

#include <L0/Print.h>
#include <L0/Ahmidas.h>
#include <L0/Core/Propagator.h>
#include <L0/Core/Correlator.h>
#include <L1/Tool/IO.h>
#include <L2/Contract/Meson.h>
#include <L2/Input/FileReader.h>

using namespace std;

int main(int narg, char **arg)
{
  if(narg<2)
    {
      cerr<<"Use "<<arg[0]<<" inputfile"<<endl;
      cerr<<"with:"<<endl;
      cerr<<" L,T"<<endl;
      cerr<<" Kappa"<<endl;
      cerr<<" Tslice of source"<<endl;
      cerr<<" Norm of the source (for normalization)"<<endl;
      cerr<<" Gauge configuration file"<<endl;
      cerr<<" Number of contractions"<<endl;
      cerr<<" List of contraction pairs"<<endl;
      cerr<<" Number of propagators"<<endl;
      cerr<<" Base name of the propagator, mass, theta"<<endl;
      cerr<<" Number of combinations"<<endl;
      cerr<<" Prop1, if1, prop2, if2"<<endl;
      exit(1);
    }

  ifstream ifile(arg[1]);

  int L,T;
  ifile>>L>>T;
  Base::Weave weave(L,T);
  bool isr=weave.isRoot();

  double kappa;
  ifile>>kappa;
  if(isr) cout<<"kappa: "<<kappa<<endl;
  int tss;
  ifile>>tss;
  if(isr) cout<<"tss: "<<tss<<endl;
  double norm;
  ifile>>norm;
  if(isr) cout<<"norm: "<<norm<<endl;

  //Load gauge field
  string gauge_file;
  Core::Field<QCD::Gauge> gauge_field(L,T);
  ifile>>gauge_file;
  Tool::IO::load(&gauge_field,gauge_file,Tool::IO::fileILDG);
  if(isr) cout<<"Gauge "<<gauge_file<<" loaded"<<endl;

  int ncontr;
  ifile>>ncontr;
  if(isr) cout<<"Ncontr: "<<ncontr<<endl;
  string *op1=new string[ncontr],*op2=new string[ncontr];
  std::vector<std::pair<Base::Operator,Base::Operator> > op_comb;
  for(int icontr=0;icontr<ncontr;icontr++)
    {
      ifile>>op1[icontr]>>op2[icontr];
      op_comb.push_back(std::make_pair(Tool::convertIntToOperator(atoi(op1[icontr].c_str())),Tool::convertIntToOperator(atoi(op2[icontr].c_str()))));
      if(isr) cout<<"Contraction "<<op1[icontr]<<" "<<op2[icontr]<<" added"<<endl;
    }

  int nprop;
  ifile>>nprop;
  if(isr) cout<<"Nprop to load: "<<nprop<<endl;
  const double thet=1;
  double *mass=new double[nprop];
  double *thes=new double[nprop];
  PropagatorType ***prop=new PropagatorType**[nprop];
  PropagatorType temp_prop(L,T);
  for(int iprop=0;iprop<nprop;iprop++)
    {
      string path_beg,path_end="";
      ifile>>path_beg>>mass[iprop]>>thes[iprop];
      if(isr) cout<<"file "<<iprop<<" "<<path_beg<<" "<<mass[iprop]<<" "<<thes[iprop]<<endl;

      char ind[3];
      vector<string> prop_path;
      prop[iprop]=new PropagatorType*[2];
      prop[iprop][0]=new PropagatorType(L,T);
      prop[iprop][1]=new PropagatorType(L,T);

      prop_path.clear();
      for(int ifi=0;ifi<FilesPerProp;ifi++)
	{
	  sprintf(ind,"%02d",ifi);
	  prop_path.push_back(path_beg+ind+path_end);
	}
      
      Tool::IO::load(&temp_prop,prop_path,Tool::IO::fileSCIDAC);
      double the=thes[iprop];
      temp_prop.reconstruct_doublet(*(prop[iprop][0]),*(prop[iprop][1]),gauge_field,kappa,mass[iprop],thet,the,the,the);
      
      prop[iprop][0]->rotateToPhysicalBasis(0);
      prop[iprop][1]->rotateToPhysicalBasis(1);
    }

  std::vector<Core::Correlator<Dirac::Matrix> > C2;

  int ncombo;
  ifile>>ncombo;
  ofstream fout;
  if(isr) fout.open("2pts_correlations");
  fout.precision(12);

  if(isr) cout<<"Ncombo: "<<ncombo<<endl;

  for(int icombo=0;icombo<ncombo;icombo++)
    {
      int iprop1,iprop2,if1,if2;
      ifile>>iprop1>>if1>>iprop2>>if2;
      if(isr)
	{
	  cout<<iprop1<<" "<<if1<<" "<<iprop2<<" "<<if2<<" "<<endl;
	  fout<<" # m1 = "<<mass[iprop1]<<" , r1 = "<<if1<<" , th1 = "<<thes[iprop1]
	      <<" , m2 = "<<mass[iprop2]<<" , r2 = "<<if2<<" , th2 = "<<thes[iprop2]<<endl;
	  fout<<endl;
	}

#ifdef UltraStocCase
      C2=Contract::light_meson_twopoint_ultrastochastic(*(prop[iprop1][if1]),*(prop[iprop2][if2]),op_comb);
#elif defined StocCase
      C2=Contract::light_meson_twopoint_stochastic(*(prop[iprop1][if1]),*(prop[iprop2][if2]),op_comb);
#else
      C2=Contract::light_meson_twopoint(*(prop[iprop1][if1]),*(prop[iprop2][if2]),op_comb);
#endif
      if(isr)
	{
	  for(int icontr=0;icontr<ncontr;icontr++)
	    {
	      fout<<" # op1 = "<<op_comb[icontr].first<<" , op2 = "<<op_comb[icontr].second<<endl;
	      C2[icontr].setOffset(tss);
	      C2[icontr]*=1/norm;
	      for(int t=0;t<T;t++)
		{
		  complex<double> corr=C2[icontr][t].trace();
		  fout<<corr.real()<<" "<<corr.imag()<<endl;
		}
	      fout<<endl;
	    }
	  fout<<endl;
	}
    }
  
  if(isr) fout.close();

  return EXIT_SUCCESS;
}
