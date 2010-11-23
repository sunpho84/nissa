#include <mpi.h>
#include <fstream>
#include <lemon.h>

#include "appretto.h"

using namespace std;

int main(int narg,char **arg)
{
  //basic initialization
  init_appretto();

  open_input("input");

  read_int("L",glb_size[1]);
  read_int("T",glb_size[0]);

  close_input();

  //Init the MPI grid 
  init_grid();
  
  ///////////////////////////////////////////

  spincolor *spinore=new spincolor[loc_vol];

  char filename[1024]="akatabum.00";
  read_spincolor(filename,spinore);

  //Print the spincolor
  for(int ivol=0;ivol<loc_vol;ivol++)
    for(int id1=0;id1<4;id1++)
      for(int ic1=0;ic1<3;ic1++)
        for(int im=0;im<2;im++)
          cout<<rank<<" "<<ivol<<" "<<id1<<" "<<ic1<<" "<<im<<" "<<spinore[ivol][id1][ic1][im]<<endl;

  
  ///////////////////////////////////////////

  close_appretto();

  return 0;
}
