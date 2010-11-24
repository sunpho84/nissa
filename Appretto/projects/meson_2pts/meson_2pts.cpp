#include <mpi.h>
#include <fstream>
#include <lemon.h>

#include "appretto.h"

using namespace std;

int main(int narg,char **arg)
{
  char base_filename1[1024],base_filename2[1024];

  //basic mpi initialization
  init_appretto();

  if(narg<2)
    {
      if(rank==0) cerr<<"Use: "<<arg[0]<<" input_file"<<endl;
      MPI_Abort(MPI_COMM_WORLD,1);
    }

  open_input(arg[1]);

  read_int("L",glb_size[1]);
  read_int("T",glb_size[0]);
  read_str("BaseFilename1",base_filename1);
  read_str("BaseFilename2",base_filename2);

  close_input();

  //Init the MPI grid 
  init_grid();
  
  ///////////////////////////////////////////

  colorspinspin *spinore1=new colorspinspin[loc_vol];
  colorspinspin *spinore2=new colorspinspin[loc_vol];

  complex *contr=new complex[glb_size[0]];

  read_colorspinspin(base_filename1,spinore1);
  read_colorspinspin(base_filename2,spinore2);
  
  //take initial time
  double tic;
  if(debug)
    {
      MPI_Barrier(cart_comm);
      tic=MPI_Wtime();
    }
  contract_two_spinspins(contr,base_gamma[5],spinore1,base_gamma[5],spinore2,-1,-1);
  //take final time
  double tac;
  if(debug)
    {
      MPI_Barrier(cart_comm);
      tac=MPI_Wtime();
      if(rank==0) cout<<"Time elapsed in contracting: "<<tac-tic<<" s"<<endl;
    }
  
  int spat_vol=glb_size[1]*glb_size[1]*glb_size[1];
  if(rank==0) for(int t=0;t<glb_size[0];t++) cout<<t<<" "<<contr[t][0]/spat_vol<<" "<<contr[t][1]/spat_vol<<endl;

  delete[] spinore1;
  delete[] spinore2;

  ///////////////////////////////////////////

  close_appretto();

  return 0;
}
