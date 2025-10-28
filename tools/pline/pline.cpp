#include "nissa.h"

//This works in serial

void Pline_serial(su3 *Pline, quad_su3 *conf)
{
  su3 *U0=(su3*)malloc(sizeof(su3)*(loc_size[0]+1));
  su3 *plf=(su3*)malloc(sizeof(su3)*(loc_size[0]+1));
  su3 *plb=(su3*)malloc(sizeof(su3)*(loc_size[0]+1));

  int iX[4];
  int L=loc_size[1];
  int T=loc_size[0];
  int X;

  for (iX[1]=0; iX[1]<L; iX[1]++)
    for (iX[2]=0; iX[2]<L; iX[2]++)
      for (iX[3]=0; iX[3]<L; iX[3]++)
	{
	  for (iX[0]=0; iX[0]<T; iX[0]++)
	    {
	      X=loclxOfCoord(iX);
	      unsafe_su3_hermitian(U0[iX[0]+1],conf[X][0]);
	    }
	  su3_copy(plf[1],U0[1]);
	  for (int t=2; t<=T; t++) su3_prod_su3(plf[t],U0[t],plf[t-1]);
	  for (iX[0]=0; iX[0]<T; iX[0]++)
	    {
	      X=loclxOfCoord(iX);
	      su3_copy(U0[iX[0]+1],conf[X][0]);
	    }
	  su3_copy(plb[T],U0[T]);
	  for (int t=1; t<T; t++)
	    {
	      int tback=T-t;
	      su3_prod_su3(plb[tback],U0[tback],plb[tback+1]);
	    }
	  iX[0]=0;
	  X=loclxOfCoord(iX);
	  su3_put_to_id(Pline[X]);
	  for (iX[0]=1; iX[0]<T; iX[0]++)
	    {
	      int t=iX[0];
	      X=loclxOfCoord(iX);
	      if( t <= ( T/2 - 1 ) ) su3_copy(Pline[X],  plf[ t ]);
	      else  su3_copy(Pline[X],plb[t+1]);
	    }	      
	}
}

void Pline_serial_forward(su3 *Pline, quad_su3 *conf)
{
  su3 *U0=(su3*)malloc(sizeof(su3)*(loc_size[0]+1));
  su3 *plf=(su3*)malloc(sizeof(su3)*(loc_size[0]+1));

  int iX[4];
  int L=loc_size[1];
  int T=loc_size[0];
  int X;

  for (iX[1]=0; iX[1]<L; iX[1]++)
    for (iX[2]=0; iX[2]<L; iX[2]++)
      for (iX[3]=0; iX[3]<L; iX[3]++)
	{
	  for (iX[0]=0; iX[0]<T; iX[0]++)
	    {
	      X=loclxOfCoord(iX);
	      unsafe_su3_hermitian(U0[iX[0]+1],conf[X][0]);
	    }
	  su3_copy(plf[1],U0[1]);
	  for (int t=2; t<=T; t++) su3_prod_su3(plf[t],U0[t],plf[t-1]);
	  
	  iX[0]=0;
	  X=loclxOfCoord(iX);
	  su3_put_to_id(Pline[X]);
	  for (iX[0]=1; iX[0]<T; iX[0]++)
	    {
	      int t=iX[0];
	      X=loclxOfCoord(iX);
	      su3_copy(Pline[X],  plf[ t ]);
	    }
	}
}

void Pline_serial_backward(su3 *Pline, quad_su3 *conf)
{
  su3 *U0=(su3*)malloc(sizeof(su3)*(loc_size[0]+1));
  su3 *plb=(su3*)malloc(sizeof(su3)*(loc_size[0]+1));

  int iX[4];
  int L=loc_size[1];
  int T=loc_size[0];
  int X;
   
  for (iX[1]=0; iX[1]<L; iX[1]++)
    for (iX[2]=0; iX[2]<L; iX[2]++)
      for (iX[3]=0; iX[3]<L; iX[3]++)
	{
	  for (iX[0]=0; iX[0]<T; iX[0]++)
	    {
	      X=loclxOfCoord(iX);
	      su3_copy(U0[iX[0]+1],conf[X][0]);
	    }
	  su3_copy(plb[T],U0[T]);
	  for (int t=1; t<T; t++)
	    {
	      int tback=T-t;
	      su3_prod_su3(plb[tback],U0[tback],plb[tback+1]);
	    }
	  for (iX[0]=T-1; iX[0]>=0; iX[0]--)
	    {
	      int t=iX[0];
	      X=loclxOfCoord(iX);
	      su3_copy(Pline[X],plb[t+1]);
	    }
	}
}

int main(int narg,char **arg)
{
  char filename[1024];

  //basic mpi initialization
  initNissa();

  if(narg<2 && rank==0)
    {
      fprintf(stderr,"Use: %s input_file\n",arg[0]);
      fflush(stderr);
      MPI_Abort(MPI_COMM_WORLD,1);
    }

  open_input(arg[1]);

  int L,T;
  read_str_int("L",&L);
  read_str_int("T",&T);
  read_str_str("Filename",filename,1024);

  close_input();

  //Init the MPI grid 
  initGrid(T,L);

  ///////////////////////////////////////////

  quad_su3 *conf=(quad_su3*)malloc(sizeof(quad_su3)*(loc_vol));
  read_ildg_gauge_conf(conf,filename);

  su3 *Pl=(su3*)malloc(sizeof(su3)*(loc_vol));
  //Pline_serial_forward(Pl, conf);
  //Pline_forward(Pl, conf);
  //Pline_serial(Pl, conf);
  //Pline(Pl, conf);
  //Pline_serial_backward(Pl,conf);
  Pline_backward(Pl,conf);
  
  int ic_in, ic_out; //auxiliar color indices

  for (int rank_x0=0; rank_x0<nrank_dir[0]; rank_x0++) for (int rank_x1=0; rank_x1<nrank_dir[1]; rank_x1++) for (int rank_x2=0; rank_x2<nrank_dir[2]; rank_x2++) for (int rank_x3=0; rank_x3<nrank_dir[3]; rank_x3++){
    if(rank_coord[0]==rank_x0 && rank_coord[1]==rank_x1 && rank_coord[2]==rank_x2 && rank_coord[3]==rank_x3){
  //for (int irank=0; irank<nissa_nranks ; irank++){ if(rank==irank){
	for(int ivol=0;ivol<loc_vol;ivol++){
  	printf("(t,x,y,z)=(%d,%d,%d,%d)\n",glb_coord_of_loclx[ivol][0],glb_coord_of_loclx[ivol][1],glb_coord_of_loclx[ivol][2],glb_coord_of_loclx[loc_site][3]);
	  for (ic_in=0; ic_in<3; ic_in++){
		printf("[");
		for(ic_out=0; ic_out<3; ic_out++) printf("%8.8f +I %8.8f \t",Pl[ivol][ic_in][ic_out][0],Pl[ivol][ic_in][ic_out][1]);
		printf("]\n");
           }
         }
	}
        MPI_Barrier(MPI_COMM_WORLD);
   }




  free(conf);
  free(Pl);
  
  ///////////////////////////////////////////

  closeNissa();

  return 0;
}
