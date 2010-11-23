#include <mpi.h>
#include <fstream>
#include <lemon.h>

#include "appretto.h"

using namespace std;

int main(int narg,char **arg)
{
  char base_filename[1024];

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
  read_str("BaseFilename",base_filename);

  close_input();

  //Init the MPI grid 
  init_grid();
  
  ///////////////////////////////////////////

  colorspinspin *spinore=new colorspinspin[loc_vol];

  read_colorspinspin(base_filename,spinore);

  if(rank==0)
    for(int loc_site=0;loc_site<loc_vol;loc_site++)
      for(int ic=0;ic<3;ic++)
	{
	  for(int id_sink=0;id_sink<4;id_sink++)
	    {
	      for(int id_source=0;id_source<4;id_source++)
		cout<<spinore[loc_site][ic][id_source][id_sink][0]<<","<<spinore[loc_site][ic][id_source][id_sink][0]<<"\t";
	      cout<<endl;
	    }
	  cout<<endl;
	}
	      
  delete[] spinore;

  ///////////////////////////////////////////

  close_appretto();

  return 0;
}
