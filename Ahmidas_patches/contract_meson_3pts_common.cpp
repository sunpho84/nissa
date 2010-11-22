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
      cerr<<" Flavor of the spectator (0,1)"<<endl;
      cerr<<" Kappa"<<endl;
      cerr<<" Tslice of source"<<endl;
      cerr<<" Norm of the source (for normalization)"<<endl;
      cerr<<" Gauge configuration file"<<endl;
      cerr<<" Number of contractions"<<endl;
      cerr<<" List of contraction pairs for: source, insertion"<<endl;
      cerr<<" Number of S0 propagators"<<endl;
      cerr<<" Base name of the propagator, mass, theta"<<endl;
      cerr<<" Number of S1 propagators"<<endl;
      cerr<<" Base name of the propagator, mass, theta"<<endl;
      cerr<<" Number of combinations between S0 S1"<<endl;
      cerr<<" Prop0, prop1"<<endl;
      cerr<<" Original source"<<endl;
      exit(1);
    }

  ifstream ifile(arg[1]);

  int L,T;
  ifile>>L>>T;
  Base::Weave weave(L,T);
  bool isr=weave.isRoot();
  Dirac::Gamma<5> gamma5;

  double kappa;
  ifile>>kappa;
  if(isr) cout<<"kappa: "<<kappa<<endl;
  bool fS0;
  ifile>>fS0;
  if(isr) cout<<"flavor of S0: "<<fS0<<endl;
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
  string *op0=new string[ncontr],*op1=new string[ncontr];
  std::vector<std::pair<Base::Operator,Base::Operator> > op_comb;
  for(int icontr=0;icontr<ncontr;icontr++)
    {
      ifile>>op0[icontr]>>op1[icontr];
      op_comb.push_back(std::make_pair(Tool::convertIntToOperator(atoi(op0[icontr].c_str())),Tool::convertIntToOperator(atoi(op1[icontr].c_str()))));
      if(isr) cout<<"Contraction "<<op0[icontr]<<" "<<op1[icontr]<<" added"<<endl;
    }

  PropagatorType temp_prop(L,T);
  const double thet=1;

  int npropS0;
  ifile>>npropS0;
  if(isr) cout<<"Nprop S0 to load: "<<npropS0<<endl;
  double *massS0=new double[npropS0];
  double *thesS0=new double[npropS0];
  PropagatorType **propS0=new PropagatorType*[npropS0];
  for(int ipropS0=0;ipropS0<npropS0;ipropS0++)
    {
      string path_beg,path_end="";
      ifile>>path_beg>>massS0[ipropS0]>>thesS0[ipropS0];
      if(isr) cout<<"file "<<ipropS0<<" "<<path_beg<<" "<<massS0[ipropS0]<<" "<<thesS0[ipropS0]<<endl;

      char ind[3];
      vector<string> prop_path;
      propS0[ipropS0]=new PropagatorType(L,T);

      prop_path.clear();
      for(int ifi=0;ifi<FilesPerProp;ifi++)
	{
	  sprintf(ind,"%02d",ifi);
	  prop_path.push_back(path_beg+ind+path_end);
	}
      
      Tool::IO::load(&temp_prop,prop_path,Tool::IO::fileSCIDAC);
      double the=thesS0[ipropS0];
      if(fS0==0) (*(propS0[ipropS0]))=temp_prop.applyDiracOperator(gauge_field,kappa,+massS0[ipropS0],thet,the,the,the);
      else       (*(propS0[ipropS0]))=temp_prop.applyDiracOperator(gauge_field,kappa,-massS0[ipropS0],thet,the,the,the);
      (*(propS0[ipropS0])).rightMultiply(gamma5); //this is required because D=Q*g5 and g5 is not put by applyDiracOperator
      
      propS0[ipropS0]->rotateToPhysicalBasis(fS0);
    }

  int npropS1;
  ifile>>npropS1;
  if(isr) cout<<"Nprop S1 to load: "<<npropS1<<endl;
  double *massS1=new double[npropS1];
  double *thesS1=new double[npropS1];
  PropagatorType **propS1=new PropagatorType*[npropS1];
  for(int ipropS1=0;ipropS1<npropS1;ipropS1++)
    {
      string path_beg,path_end="";
      ifile>>path_beg>>massS1[ipropS1]>>thesS1[ipropS1];
      if(isr) cout<<"file "<<ipropS1<<" "<<path_beg<<" "<<massS1[ipropS1]<<" "<<thesS1[ipropS1]<<endl;

      char ind[3];
      vector<string> prop_path;
      propS1[ipropS1]=new PropagatorType(L,T);

      prop_path.clear();
      for(int ifi=0;ifi<FilesPerProp;ifi++)
	{
	  sprintf(ind,"%02d",ifi);
	  prop_path.push_back(path_beg+ind+path_end);
	}
      
      Tool::IO::load(&temp_prop,prop_path,Tool::IO::fileSCIDAC);
      double the=thesS1[ipropS1];
      //Now apply the dirac operator. It has to be the opposite of what
      //you want in the s1 line, so the same of s0_flav (0->-,1->+)
      if(fS0==0) (*(propS1[ipropS1]))=temp_prop.applyDiracOperator(gauge_field,kappa,-massS1[ipropS1],thet,the,the,the);
      else       (*(propS1[ipropS1]))=temp_prop.applyDiracOperator(gauge_field,kappa,+massS1[ipropS1],thet,the,the,the);
      (*(propS1[ipropS1])).rightMultiply(gamma5); //this is required because D=Q*g5 and g5 is not put by applyDiracOperator

      //Now rotate on the source as s0 (0->-,1->+)                                                                                
      PropagatorType temp(L,T);
#ifndef UltraStocCase
      temp=(*(propS1[ipropS1]));
      temp.isolate();
      if(fS0==0) temp*=std::complex<double>(0,-1);
      else       temp*=std::complex<double>(0,+1);
      temp*=(gamma5);
      (*(propS1[ipropS1]))+=temp;
#endif

      //Now rotate on the sink as s1 (0->+,1->-)
      temp=(*(propS1[ipropS1]));
      temp.isolate();
      if(fS0==0) temp*=std::complex<double>(0,+1);
      else       temp*=std::complex<double>(0,-1);
      temp.rightMultiply(gamma5);
      (*(propS1[ipropS1]))+=temp;

      //Multiply by 1/sqrt(2)**2                                                                                                  
#ifndef UltraStocCase
      (*(propS1[ipropS1]))*=0.5;
#else
      (*(propS1[ipropS1]))*=1/sqrt(2);
#endif
    }

  std::vector<Core::Correlator<Dirac::Matrix> > C2;

  ofstream fout;
  fout.precision(12);
  if(isr) fout.open("3pts_correlations");

  int ncombo;
  ifile>>ncombo;
  if(isr) cout<<"Ncombo: "<<ncombo<<endl;

  for(int icombo=0;icombo<ncombo;icombo++)
    {
      int iprop0,iprop1;
      ifile>>iprop0>>iprop1;
      if(isr)
        {
          cout<<iprop0<<" "<<iprop1<<" "<<endl;
          fout<<" # r_spec = "<<fS0<<" , m0 = "<<massS0[iprop0]<<" , th0 = "<<thesS0[iprop0]
	                           <<" , mSeq = "<<massS1[iprop1]<<" , thSeq = "<<thesS1[iprop1]<<endl;
        }
#ifdef UltraStocCase
      C2=Contract::light_meson_twopoint_ultrastochastic(*(propS1[iprop1]),*(propS0[iprop0]),op_comb);
#elif defined StocCase
      C2=Contract::light_meson_twopoint_stochastic(*(propS1[iprop1]),*(propS0[iprop0]),op_comb);
#else
      C2=Contract::light_meson_twopoint(*(propS1[iprop1]),*(propS0[iprop0]),op_comb);
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

  std::vector<std::pair<Base::Operator,Base::Operator> > op_2pts_p5p5;
  op_2pts_p5p5.push_back(std::make_pair(Tool::convertIntToOperator(-1),Tool::convertIntToOperator(5)));

  vector<string> prop_path;
  string path_source,path_2pts_check;
  ifile>>path_source>>path_2pts_check;
  prop_path.clear();
  char ind[3];
  for(int ifi=0;ifi<FilesPerProp;ifi++)
    {
      sprintf(ind,"%02d",ifi);
      prop_path.push_back(path_source+ind);
    }
      
  Tool::IO::load(&temp_prop,prop_path,Tool::IO::fileSCIDAC);

  if(isr)
    {
      fout.clear();
      fout.open("2pts_check");
      fout.precision(12);
    }

  for(int ipropS1=0;ipropS1<npropS1;ipropS1++)
    {
#ifdef UltraStocCase
      C2=Contract::light_meson_twopoint_ultrastochastic(*(propS1[ipropS1]),temp_prop,op_2pts_p5p5);
#elif defined StocCase
      C2=Contract::light_meson_twopoint_stochastic(*(propS1[ipropS1]),temp_prop,op_2pts_p5p5);
#else
      C2=Contract::light_meson_twopoint(*(propS1[ipropS1]),temp_prop,op_2pts_p5p5);
#endif
      if(isr)
	{
	  fout<<" # m1 = "<<massS1[ipropS1]<<" , th1 = "<<thesS1[ipropS1]<<endl;
	  
	  C2[0].setOffset(tss);
	  C2[0]*=1/norm;
	  for(int t=0;t<T;t++)
	    {
	      complex<double> corr=C2[0][t].trace();
	      fout<<corr.real()<<" "<<corr.imag()<<endl;
	    }
	  fout<<endl;
	}
    }

  fout.close();
  
  return EXIT_SUCCESS;
}
