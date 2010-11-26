#include <mpi.h>
#include <fstream>
#include <lemon.h>

#include "appretto.h"

using namespace std;

/*
             SEQ	                  	 
            .....                        .....           
          ..     ..                    ..     ..         
         .         .                  .         .        
    op  X           X  g5   =    op  X           .  S1 
         .         .                  .         .        
          ..     ..                    ..     ..         
            .....                        .....           
              S0                                    
*/

//This function contract a source with a sequential spinor putting the passed list of operators
void contract_with_source(complex *corr,colorspinspin *S1,int *list_op,colorspinspin *source,int ncontr,int fprop,int rprop)
{
  //Temporary vector for the internal matrices
  dirac_matr t1[ncontr],t2[ncontr];

  for(int icontr=0;icontr<ncontr;icontr++)
    {
      
      //Init the second 
      t2[icontr]=base_gamma[0];
      
      //Remind that D- rotates as 1+ig5, but D-^-1 rotates as 1-ig5,

      //f1 < 1: do rotation for a quark propagator
      //f1 = 1: do nothing (physical basis already)
      //f1 > 1: do rotation for a (charged) sequential propagator

      if(fprop<1)
	switch(rprop)
	  {
	  case 0: //This is D-^-1
	    dirac_prod(t1[icontr], t1[icontr],Pminus);
	    dirac_prod(t2[icontr], Pminus,t2[icontr]);
	    break;
	  case 1: //This is D+^-1
	    dirac_prod(t1[icontr], t1[icontr],Pplus);
	    dirac_prod(t2[icontr], Pplus,t2[icontr]);
	    break;
	  }

      if(fprop>1)
        switch(rprop)
          {
          case 0: //This is D-^-1
            dirac_prod(t1[icontr], t1[icontr],Pminus);
            dirac_prod(t2[icontr], Pplus,t2[icontr]);
            break;
          case 1: //This is D+^-1
            dirac_prod(t1[icontr], t1[icontr],Pplus);
            dirac_prod(t2[icontr], Pminus,t2[icontr]);
            break;
          }

    }

  //Call for the routine which does the real contraction
  trace_g_s_g_sdag(corr,t1,S1,t2,source,ncontr);
}

int main(int narg,char **arg)
{
  //Basic mpi initialization
  init_appretto();

  //take initial time
  double tic;
  if(debug)
    {
      MPI_Barrier(cart_comm);
      tic=MPI_Wtime();
    }

  if(narg<2)
    {
      if(rank==0) cerr<<"Use: "<<arg[0]<<" input_file"<<endl;
      MPI_Abort(MPI_COMM_WORLD,1);
    }

  open_input(arg[1]);

  //Read the volume
  read_int("L",glb_size[1]);
  read_int("T",glb_size[0]);

  //Init the MPI grid 
  init_grid();

  //allocate the source and the prop
  colorspinspin *source=new colorspinspin[loc_vol];
  colorspinspin *S1=new colorspinspin[loc_vol];
  if(source==NULL or S1==NULL)
    {
      if(rank==0) cerr<<"Error! Not enoug memory to allocate source and spinor!"<<endl;
      MPI_Abort(MPI_COMM_WORLD,1);
    }

  //Read the time location of the source
  int twall;
  read_int("TWall",twall);

  //Read the source
  char path[1024];
  read_str("Source",path);
  read_colorspinspin(path,source);

  //Read the number of contractions
  int ncontr;
  read_int("NContr",ncontr);
  if(rank==0) cout<<"Number of contractions: "<<ncontr<<endl;

  //Initialize the list of correlations and the list of operators
  complex contr[ncontr][glb_size[0]];
  int *op=new int[ncontr]; 
  for(int icontr=0;icontr<ncontr;icontr++)
    {
      //Read the operator pairs
      read_int(op[icontr]);

      if(rank==0 and debug) cout<<" contr."<<icontr<<" "<<op[icontr]<<endl;
    }

  //Read the number of the props
  int nprop;
  read_int("NProp",nprop);

  //Inizialization of output
  ofstream fout;
  fout.precision(16);
  fout.width(16);
  if(rank==0) fout.open("two_points_check");
  int spat_vol=glb_size[1]*glb_size[1]*glb_size[1];

  //Loop over propagators
  double mass_spec,mass1;
  double theta_spec,theta1;
  int phys;
  int r1;
  for(int iprop=0;iprop<nprop;iprop++)
    {

      //Read the name, mass, theta and r for the list of the propagators
      read_str(path,1024);
      read_double(mass_spec);
      read_double(theta_spec);
      read_double(mass1);
      read_double(theta1);
      read_int(phys);
      read_int(r1);

      if(debug and rank==0)
	cout<<" prop."<<iprop<<" "<<path<<", m_spec="<<mass_spec<<" th_spec="<<theta_spec<<", m1="<<mass1<<" th1="<<theta1<<" phys="<<phys<<" r1="<<r1<<endl;
      
      //Read the propagator
      read_colorspinspin(path,S1);
      
      if(rank==0)
	fout<<noshowpos<<" #"
	    <<" m_spec="<<mass_spec<<" th_spec="<<theta_spec
	    <<" m1="<<mass1<<" th1="<<theta1<<" r1="<<r1<<endl;
      
      contract_with_source((complex*)contr,S1,op,source,ncontr,phys,r1);

      if(rank==0)
	{
	  fout<<endl;
	  
	  for(int icontr=0;icontr<ncontr;icontr++)
	    {
	      fout<<noshowpos<<" # "<<op[icontr]<<endl;
	      fout<<endl;
	      for(int tempt=0;tempt<glb_size[0];tempt++)
		{
		  int t=tempt+twall;
		  if(t>=glb_size[0]) t-=glb_size[0];
		  
		  fout<<showpos<<contr[icontr][t][0]/spat_vol<<"\t"<<contr[icontr][t][1]/spat_vol<<endl;
		}
	      fout<<endl;
	    }
	}
      if(rank==0) fout<<endl;
    }

  close_input();

  //take final time
  double tac;
  if(debug)
    {
      MPI_Barrier(cart_comm);
      tac=MPI_Wtime();
      if(rank==0) cout<<"Time elapsed in doing "<<nprop*ncontr<<" contraction: "<<tac-tic<<" s"<<endl;
    }

  ///////////////////////////////////////////

 if(rank==0)
   {
     cout<<endl<<"Everything ok, exiting!"<<endl;
     fout.close();
   }

  delete[] S1;
  delete[] source;

  delete[] op;

  close_appretto();

  return 0;
}
