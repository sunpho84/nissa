#include <mpi.h>
#include <fstream>
#include <lemon.h>

#include "appretto.h"

using namespace std;


//Frontend to the gamma calculation
void meson_two_points(complex *corr,int *list_op1,colorspinspin *s1,int *list_op2,colorspinspin *s2,int ncontr,int r1=-1,int r2=-1)
{
  dirac_matr t1[ncontr],t2[ncontr];

  for(int icontr=0;icontr<ncontr;icontr++)
    {
      
      //Put the two gamma5 needed for the revert of the first spinor
      dirac_prod(t1[icontr], base_gamma[list_op1[icontr]],base_gamma[5]);
      dirac_prod(t2[icontr], base_gamma[5],base_gamma[list_op2[icontr]]);
      
      //Put the rotators to the physical basis these can be avoided by
      //putting a number different from 0 or 1 in f1 and f2
      
      //Remember that D- rotate as 1+ig5, but D-^-1 rotate ad 1-ig5,
      //morover (D^-1)^dagger rotate again as 1+ig5 (pweee!!!)
      
      switch(r1)
	{
	case 0: //This is (D-^-1)^dagger
	  dirac_prod(t1[icontr], t1[icontr],Pplus);
	  dirac_prod(t2[icontr], Pplus,t2[icontr]);
	  break;
	case 1:  //This is (D+^-1)^dagger
	  dirac_prod(t1[icontr], t1[icontr],Pminus);
	  dirac_prod(t2[icontr], Pminus,t2[icontr]);
	  break;
	}
      
      switch(r2)
	{
	case 0:  //This is D-^-1
	  dirac_prod(t2[icontr], t2[icontr],Pminus);
	  dirac_prod(t1[icontr], Pminus,t1[icontr]);
	  break;
	case 1:  //This is D+^-1
	  dirac_prod(t2[icontr], t2[icontr],Pplus);
	  dirac_prod(t1[icontr], Pplus,t1[icontr]);
	  break;
	}
    }
  
  trace_g_sdag_g_s(corr,t1,s1,t2,s2,ncontr);
}

int main(int narg,char **arg)
{
  //basic mpi initialization
  init_appretto();

  if(narg<2)
    {
      if(rank==0) cerr<<"Use: "<<arg[0]<<" input_file"<<endl;
      MPI_Abort(MPI_COMM_WORLD,1);
    }

  open_input(arg[1]);

  //Read the volume
  read_int("L",glb_size[1]);
  read_int("T",glb_size[0]);

  //Read the two spinors
  char base_filename1[1024];
  char base_filename2[1024];
  read_str("BaseFilename1",base_filename1);
  read_str("BaseFilename2",base_filename2);

  //Read the number of contractions
  int ncontr;
  read_int("NContr",ncontr);

  //Initialize the list of correlations and the list of operators
  complex contr[ncontr][glb_size[0]];
  int *op1=new int[ncontr]; 
  int *op2=new int[ncontr]; 
  for(int icontr=0;icontr<ncontr;icontr++)
    {
      if(rank==0) cout<<"addr:"<<contr<<endl;
      //read the operators pairs
      read_int(op1[icontr]);
      read_int(op2[icontr]);

      if(rank==0) cout<<" contr. "<<icontr<<" "<<op1[icontr]<<" "<<op2[icontr]<<endl;
    }

  close_input();

  //Init the MPI grid 
  init_grid();
  
  ///////////////////////////////////////////

  //read the spinors
  colorspinspin *spinor1=new colorspinspin[loc_vol];
  colorspinspin *spinor2=new colorspinspin[loc_vol];
  read_colorspinspin(base_filename1,spinor1);
  read_colorspinspin(base_filename2,spinor2);
  
  //take initial time
  double tic;
  if(debug)
    {
      MPI_Barrier(cart_comm);
      tic=MPI_Wtime();
    }
  meson_two_points((complex*)contr,op1,spinor1,op2,spinor2,ncontr);
  //take final time
  double tac;
  if(debug)
    {
      MPI_Barrier(cart_comm);
      tac=MPI_Wtime();
      if(rank==0) cout<<"Time elapsed in contracting: "<<tac-tic<<" s"<<endl;
    }
  
  int spat_vol=glb_size[1]*glb_size[1]*glb_size[1];
  if(rank==0)
    for(int icontr=0;icontr<ncontr;icontr++)
      for(int t=0;t<glb_size[0];t++)
	cout<<t<<" "<<contr[icontr][t][0]/spat_vol<<" "<<contr[icontr][t][1]/spat_vol<<endl;

  ///////////////////////////////////////////

  delete[] spinor1;
  delete[] spinor2;

  delete[] op1;
  delete[] op2;

  close_appretto();

  return 0;
}
