#include "nissa.h"

int main(int narg,char **arg)
{
  char filename[1024];

  //basic mpi initialization
  init_nissa();

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
  //Init the MPI grid 
  init_grid(T,L);
  
  read_str_str("Filename",filename,1024);

  double theta[4];
  read_str_double("ThetaTXYZ",&(theta[0]));
  read_double(&(theta[1]));
  read_double(&(theta[2]));
  read_double(&(theta[3]));

  if(rank==0)
    printf("Boundary condition (in multiples of 2Pi : %f %f %f %f\n",theta[0],theta[1],theta[2],theta[3]);

  int Nindexmom;
  read_str_int("NindexMom",&Nindexmom);
  if(rank==0) printf("Nindexmom = %d\n",Nindexmom);

  int ***iP=(int***)malloc(sizeof(int*)*(Nindexmom));
  for(int imom=0;imom<Nindexmom;imom++) iP[imom]=(int**)malloc(sizeof(int*)*(4));
  for(int imom=0;imom<Nindexmom;imom++) for (int idir=0; idir<4; idir++) iP[imom][idir]=(int*)malloc(sizeof(int*)*(2));

  for(int imom=0;imom<Nindexmom;imom++){
   read_str_int("iPT",&(iP[imom][0][0]));
   read_int(&(iP[imom][0][1]));
   read_str_int("iPX",&(iP[imom][1][0]));
   read_int(&(iP[imom][1][1]));
   read_str_int("iPY",&(iP[imom][2][0]));
   read_int(&(iP[imom][2][1]));
   read_str_int("iPZ",&(iP[imom][3][0]));
   read_int(&(iP[imom][3][1]));
  }
  if(rank==0){
   for(int imom=0;imom<Nindexmom;imom++){
        printf("Moment set %d\n",imom);
	printf("iPT_i =%d iPT_f=%d ; iPX_i =%d iPX_f=%d ; iPY_i =%d iPY_f=%d ; iPZ_i =%d iPZ_f=%d\n", iP[imom][0][0],iP[imom][0][1],iP[imom][1][0],iP[imom][1][1],iP[imom][2][0],iP[imom][2][1],iP[imom][3][0],iP[imom][3][1]);
    }
  }


  close_input();



  int Nmom_per_set[Nindexmom];
  //Comprueba si realmente necesito tener esta información en todos los nodos
  for (int irank=0; irank<nissa_nranks; irank++){
    if(rank==irank)   for(int imom=0;imom<Nindexmom;imom++){
		Nmom_per_set[imom]=1;
		for(int idir=0; idir<4; idir++) Nmom_per_set[imom]=Nmom_per_set[imom]*(iP[imom][idir][1]-iP[imom][idir][0]+1);
     }
  MPI_Barrier(MPI_COMM_WORLD);
  }

  if (rank==0) {
	for (int indexmom=0; indexmom<Nindexmom; indexmom++) printf("Nmom in set  %d : %d\n", indexmom, Nmom_per_set[indexmom]);
  }
  
  ///////////////////////////////////////////
  spincolor *spinore=(spincolor*)malloc(sizeof(spincolor)*loc_vol);

  read_spincolor(spinore,filename);

  spincolor **spinoreFT=(spincolor**)malloc(sizeof(spincolor*)*(Nindexmom)); 
  for(int indexmom=0;indexmom<Nindexmom;indexmom++) spinoreFT[indexmom]=(spincolor*)malloc(sizeof(spincolor)*(Nmom_per_set[indexmom]));
 
  for(int indexmom=0;indexmom<Nindexmom;indexmom++)	spincolor_FT(spinore,spinoreFT[indexmom],theta,iP[indexmom],Nmom_per_set[indexmom]);

      double **P2=(double**)malloc(sizeof(double)*(Nindexmom));
      for(int indexmom=0;indexmom<Nindexmom;indexmom++) P2[indexmom]=(double*)malloc(sizeof(double)*(Nmom_per_set[indexmom]));
      double **SinP2=(double**)malloc(sizeof(double)*(Nindexmom));
      for(int indexmom=0;indexmom<Nindexmom;indexmom++) SinP2[indexmom]=(double*)malloc(sizeof(double)*(Nmom_per_set[indexmom]));
      double ***P=(double***)malloc(sizeof(double)*(Nindexmom));
      for(int indexmom=0;indexmom<Nindexmom;indexmom++) P[indexmom]=(double**)malloc(sizeof(double)*(4));
      for(int indexmom=0;indexmom<Nindexmom;indexmom++) for (int idir=0 ; idir<4; idir++) P[indexmom][idir]=(double*)malloc(sizeof(double)*(Nmom_per_set[indexmom]));
      double ***SinP=(double***)malloc(sizeof(double)*(Nindexmom));
      for(int indexmom=0;indexmom<Nindexmom;indexmom++) SinP[indexmom]=(double**)malloc(sizeof(double)*(4));
      for(int indexmom=0;indexmom<Nindexmom;indexmom++) for (int idir=0 ; idir<4; idir++) SinP[indexmom][idir]=(double*)malloc(sizeof(double)*(Nmom_per_set[indexmom]));
      double **SinP4=(double**)malloc(sizeof(double)*(Nindexmom));
      for(int indexmom=0;indexmom<Nindexmom;indexmom++) SinP4[indexmom]=(double*)malloc(sizeof(double)*(Nmom_per_set[indexmom]));

      if (rank==0)
      {
      	for(int indexmom=0;indexmom<Nindexmom;indexmom++)  
	{
		Momentum(iP[indexmom],theta,P2[indexmom],SinP2[indexmom],P[indexmom],SinP[indexmom],SinP4[indexmom],Nmom_per_set[indexmom]); 
       		for (int imom=0; imom<Nmom_per_set[indexmom]; imom++) 
		{
			printf("----------------------------------------------------------\n");
			printf("2piP/L=(%8.8f,%8.8f,%8.8f,%8.8f) (2piP/L)^2=%8.8f Sin(2piP/L)=(%8.8f,%8.8f,%8.8f,%8.8f) Sin2(2piP/L)=%8.8f  Sin4(2piP/L)=%8.8f\n",P[indexmom][0][imom],P[indexmom][1][imom],P[indexmom][2][imom],P[indexmom][3][imom],P2[indexmom][imom],SinP[indexmom][0][imom],SinP[indexmom][1][imom],SinP[indexmom][2][imom],SinP[indexmom][3][imom],SinP2[indexmom][imom],SinP4[indexmom][imom]);
			for (int is=0; is<4; is++) for (int ic=0; ic<3; ic++) printf("is=%d ic=%d %8.8f +I %8.8f\n",is,ic,spinoreFT[indexmom][imom][is][ic][0],spinoreFT[indexmom][imom][is][ic][1]);
		}
	}
      }

  ///////////////////////////////////////////

  close_nissa();

  return 0;
}
