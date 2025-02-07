#include "appretto.h"

int main(int narg,char **arg)
{
  //basic mpi initialization
  init_appretto();

  if(narg<7 && rank==0)
    {
      fprintf(stderr,"Use: %s, L T Tslice filein fileout OUT_prec[32,64]\n",arg[0]);
      fflush(stderr);
      MPI_Abort(MPI_COMM_WORLD,1);
    }

  glb_size[1]=atoi(arg[1]);
  glb_size[0]=atoi(arg[2]);
  
  int tslice=atoi(arg[3]);

  //Init the MPI grid 
  initGrid();
  
  ///////////////////////////////////////////
  
  spincolor *spinore=(spincolor*)malloc(sizeof(spincolor)*loc_vol);
  dirac_matr gdirac=base_gamma[5];
  complex sp_temp[4][3];
  char filename[1024];
  
  //Put 0 on all timeslice different from tslice, multiply by a dirac matrix gdirac
  //and save
  for(int id1=0;id1<4;id1++)
    {
      sprintf(filename,"%s.0%d",arg[4],id1);
      read_spincolor(spinore,filename);

      NISSA_LOC_VOL_LOOP(ivol)
	{
	  if(glb_coord_of_loclx[ivol][0]!=tslice)
	    for(int id2=0;id2<4;id2++)
	      for(int ic1=0;ic1<3;ic1++)
		for(int im=0;im<2;im++)
		  spinore[ivol][id2][ic1][im]=0;
	  else
	    {
	      for(int id2=0;id2<4;id2++)
		for(int ic1=0;ic1<3;ic1++)
		  unsafe_complex_prod(sp_temp[id2][ic1],gdirac.entr[id2],spinore[ivol][gdirac.pos[id2]][ic1]);
	      for(int id2=0;id2<4;id2++)
		for(int ic1=0;ic1<3;ic1++)
		  for(int im=0;im<2;im++)
		    spinore[ivol][id2][ic1][im]=sp_temp[id2][ic1][im];
	    }
	}
      
      sprintf(filename,"%s.0%d",arg[5],id1);
      write_spincolor(filename,spinore,atoi(arg[6]));
    }
  
  free(spinore);
  
  ///////////////////////////////////////////

  close_appretto();

  return 0;
}
