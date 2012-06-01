
// This program calculates a list of three-point functions for the BK 
// contracting all the propagators in the four lists if and only if 
// the sum of the r's is odd and the masses of the first and the third lists
// as well as those of the second and the fourth lists are equal.
// The first and second lists refer to the left side (see below), while 
// the third and the fourth lists correspond to the right side 
// (different source located at T/2 away from the left one).
// The first and third lists are uploaded in blocks of doublets as large as possible.
// The second list is uploaded one by one.
// The fourth propagator is uploaded by doublets.
//
// ***** The third list MUST contain the same masses, theta's and flavors of the first list in the same order *****
//
// ***** The fourth list MUST contain the same masses, theta's and flavors of the second list in the same order *****
//
// ***** The flavors must be ordered with r=0 first *****
// 

//////////////////////////////////////////////////////////////////////////////////////////////////////
// Three-point disconnected part (Mezzotos)
//               _                                          _       
//   Sum_x{ Tr [ S(md;x,tL) gsource S(ms,x,tL) OP_L  ] Tr [ S(md;x,tR) gsource S(ms,x,tR) OP_R  ]  }      
//
// Three-point connected part (Otto)
//
//   Sum_x{ Tr [ S(md;x,tL) gsource S(ms,x,tL) OP_L  S(md;x,tR) gsource S(ms,x,tR) OP_R  ]  }    
//        _
// where  S is the revert propagator
//   _     +
//   S=g5 S g5
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
//This is the reference scheme:

/*                             
                S(ms)                           S(ms) 
                .....          			.....
              ..     ..        		      ..     ..
             .         .         (x)	     .	       .	
   gsource  X           X  OP_L        OP_R  X           X gsource
     (tL)    .         .               (tR)  .         .
              ..     ..                       ..      ..
                .....          		  	.....
                S(md)                            S(md)  

                S(ms)        S(ms)
                .....       ,.....
              ..     ..    ,.     ..   
             .         .  ,         .      
   gsource  X           X            X gsource
     (tL)    .         , .          .   (tR)
              ..     .,    ..     ..
                ....,        .....
                S(md)         S(md)

                                    
source |------>---->----->---->| sink

*/

#include <mpi.h>
#include <lemon.h>

#include "nissa.h"




//Calculate the maximum number of allocable propagators
int compute_allocable_propagators_list(int nprop_list)
{
  //Check if the nprop_list is even (flavor doublets)
  if(nprop_list != (2*(nprop_list/2))) crash("Error: number of propagators must be even");

  //Check if there is enough room for seven propagators (nprop_list = 2).
  //This is the minimal requirement for the program to be able to work.
  colorspinspin *fuf;
  fuf=(colorspinspin*)malloc(7*sizeof(colorspinspin)*loc_vol);

  if(fuf==NULL) crash("Error: not enough memory for seven propagators");
  else if(debug_lvl>1 && rank==0) printf("Ok there is enough memory to load seven propagators\n");

  free(fuf);

  //Now determine the largest number of doublets of propagators of the first and third lists (light ones) and two additional ones loadable at once.
  int nprop_max = 2 * nprop_list + 7;
  do
    {
      nprop_max-=4;
      fuf=(colorspinspin*)malloc(nprop_max*sizeof(colorspinspin)*loc_vol);
    }
  while(fuf==NULL);

  free(fuf);

  if(debug_lvl && rank==0)
    printf("Will allocate %d propagators at once \n",nprop_max);

  int nprop_max_per_list=(nprop_max-3)/2;
  return nprop_max_per_list; //return the number of allocable propagators in the first and third lists (separately)
}

//This function takes care to make the revert on the FIRST spinor, putting the needed gamma5
//It also applies the appropriate rotators to the physical basis if asked
void meson_two_points(complex **corr,int *list_op1,colorspinspin *s1,int *list_op2,colorspinspin *s2,int ncontr,int f1,int r1,int f2,int r2)
{
  //Temporary vectors for the internal gamma
  dirac_matr t1[ncontr],t2[ncontr];

  for(int icontr=0;icontr<ncontr;icontr++)
    {

      //Put the two gamma5 needed for the revert of the first spinor
      dirac_prod(&(t1[icontr]), &(base_gamma[list_op1[icontr]]),&(base_gamma[5]));
      dirac_prod(&(t2[icontr]), &(base_gamma[5]),&(base_gamma[list_op2[icontr]]));

      //Remind that D- rotates as 1+ig5, but D-^-1 rotates as 1-ig5,
      //moreover (D^-1)^dagger rotates again as 1+ig5 (pweee!!!)

      //f1 < 1: do rotation for a quark propagator
      //f1 = 1: do nothing (physical basis already)
      //f1 > 1: do rotation for a (charged) sequential propagator

      if(f1<1)
        switch(r1)
          {
          case 0: //This is (D-^-1)^dag
            dirac_prod(&(t1[icontr]), &(t1[icontr]),&Pplus);
            dirac_prod(&(t2[icontr]), &Pplus,&(t2[icontr]));
            break;
          case 1: //This is (D+^-1)^dag
            dirac_prod(&(t1[icontr]), &(t1[icontr]),&Pminus);
            dirac_prod(&(t2[icontr]), &Pminus,&(t2[icontr]));
            break;
          }

      if(f1>1)
        switch(r1)
          {
          case 0: //This is (D-^-1)^dag
            dirac_prod(&(t1[icontr]), &(t1[icontr]),&Pplus);
            dirac_prod(&(t2[icontr]), &Pminus,&(t2[icontr]));
            break;
          case 1: //This is (D+^-1)^dag
            dirac_prod(&(t1[icontr]), &(t1[icontr]),&Pminus);
            dirac_prod(&(t2[icontr]), &Pplus,&(t2[icontr]));
            break;
          }

      if(f2<1)
        switch(r2)
          {
          case 0: //This is D-^-1
            dirac_prod(&(t2[icontr]), &(t2[icontr]),&Pminus);
            dirac_prod(&(t1[icontr]), &Pminus,&(t1[icontr]));
            break;
          case 1: //This is D+^-1
            dirac_prod(&(t2[icontr]), &(t2[icontr]),&Pplus);
            dirac_prod(&(t1[icontr]), &Pplus,&(t1[icontr]));
            break;
          }

     if(f2>1)
        switch(r2)
          {
          case 0: //This is D-^-1
            dirac_prod(&(t2[icontr]), &(t2[icontr]),&Pminus);
            dirac_prod(&(t1[icontr]), &Pplus,&(t1[icontr]));
            break;
          case 1: //This is D+^-1
            dirac_prod(&(t2[icontr]), &(t2[icontr]),&Pplus);
            dirac_prod(&(t1[icontr]), &Pminus,&(t1[icontr]));
            break;
          }

    }

  //Call for the routine which does the real contraction
  trace_g_sdag_g_s(corr[0],t1,s1,t2,s2,ncontr);
}




void BK_eight_disconnected(complex **mezzo, int *list_opQ,colorspinspin *sdL,colorspinspin *ssL, colorspinspin *sdR,colorspinspin *ssR,int ncontr,int fdL,int rdL,int fsL,int rsL,int fdR,int rdR,int fsR,int rsR ){

 //Temporary vectors for the internal gamma
  dirac_matr gsource_times_g5_L[ncontr],g5_times_gQ_L[ncontr];
  dirac_matr gsource_times_g5_R[ncontr],g5_times_gQ_R[ncontr];
  int gsource=5; //Kaon-Kaon
  for(int icontr=0;icontr<ncontr;icontr++)
    {

      //Put the two gamma5 needed for the revert of the d spinor
      dirac_prod(&(gsource_times_g5_L[icontr]), &(base_gamma[gsource]),&(base_gamma[5]));
      dirac_prod(&(g5_times_gQ_L[icontr]), &(base_gamma[5]),&(base_gamma[list_opQ[icontr]]));
      dirac_prod(&(gsource_times_g5_R[icontr]), &(base_gamma[gsource]),&(base_gamma[5]));
      dirac_prod(&(g5_times_gQ_R[icontr]), &(base_gamma[5]),&(base_gamma[list_opQ[icontr]]));

     //Remind that D- rotates as 1+ig5, but D-^-1 rotates as 1-ig5,
      //moreover (D^-1)^dagger rotates again as 1+ig5 (pweee!!!)

      //f1 < 1: do rotation for a quark propagator
      //f1 = 1: do nothing (physical basis already)

      if(fdL<1)
        switch(rdL)
          {
          case 0: //This is (D-^-1)^dag
            dirac_prod(&(gsource_times_g5_L[icontr]), &(gsource_times_g5_L[icontr]),&Pplus);
            dirac_prod(&(g5_times_gQ_L[icontr]), &Pplus,&(g5_times_gQ_L[icontr]));
            break;
          case 1: //This is (D+^-1)^dag
            dirac_prod(&(gsource_times_g5_L[icontr]), &(gsource_times_g5_L[icontr]),&Pminus);
            dirac_prod(&(g5_times_gQ_L[icontr]), &Pminus,&(g5_times_gQ_L[icontr]));
            break;
          }

      if(fsL<1)
        switch(rsL)
          {
          case 0: //This is D-^-1
            dirac_prod(&(g5_times_gQ_L[icontr]), &(g5_times_gQ_L[icontr]),&Pminus);
            dirac_prod(&(gsource_times_g5_L[icontr]), &Pminus,&(gsource_times_g5_L[icontr]));
            break;
          case 1: //This is D+^-1
            dirac_prod(&(g5_times_gQ_L[icontr]), &(g5_times_gQ_L[icontr]),&Pplus);
            dirac_prod(&(gsource_times_g5_L[icontr]), &Pplus,&(gsource_times_g5_L[icontr]));
            break;
          }

      if(fdR<1)
        switch(rdR)
          {
          case 0: //This is (D-^-1)^dag
            dirac_prod(&(gsource_times_g5_R[icontr]), &(gsource_times_g5_R[icontr]),&Pplus);
            dirac_prod(&(g5_times_gQ_R[icontr]), &Pplus,&(g5_times_gQ_R[icontr]));
            break;
          case 1: //This is (D+^-1)^dag
            dirac_prod(&(gsource_times_g5_R[icontr]), &(gsource_times_g5_R[icontr]),&Pminus);
            dirac_prod(&(g5_times_gQ_R[icontr]), &Pminus,&(g5_times_gQ_R[icontr]));
            break;
          }

      if(fsR<1)
        switch(rsR)
          {
          case 0: //This is D-^-1
            dirac_prod(&(g5_times_gQ_R[icontr]), &(g5_times_gQ_R[icontr]),&Pminus);
            dirac_prod(&(gsource_times_g5_R[icontr]), &Pminus,&(gsource_times_g5_R[icontr]));
            break;
          case 1: //This is D+^-1
            dirac_prod(&(g5_times_gQ_R[icontr]), &(g5_times_gQ_R[icontr]),&Pplus);
            dirac_prod(&(gsource_times_g5_R[icontr]), &Pplus,&(gsource_times_g5_R[icontr]));
            break;
          }

    }

 //Call for the routine which does the real contraction for the Mezzotto and the otto:
  sum_trace_g_sdag_g_s_times_trace_g_sdag_g_s(mezzo,gsource_times_g5_L,sdL,g5_times_gQ_L,ssL,gsource_times_g5_R,sdR,g5_times_gQ_R,ssR,ncontr);

}




void BK_eight_connected(complex **otto, int *list_opQ,colorspinspin *sdL,colorspinspin *ssL, colorspinspin *sdR,colorspinspin *ssR,int ncontr,int fdL,int rdL,int fsL,int rsL,int fdR,int rdR,int fsR,int rsR ){

 //Temporary vectors for the internal gamma
  dirac_matr gsource_times_g5_L[ncontr],g5_times_gQ_L[ncontr];
  dirac_matr gsource_times_g5_R[ncontr],g5_times_gQ_R[ncontr];
  int gsource=5; //Kaon-Kaon
  for(int icontr=0;icontr<ncontr;icontr++)
    {

      //Put the two gamma5 needed for the revert of the d spinor
      dirac_prod(&(gsource_times_g5_L[icontr]), &(base_gamma[gsource]),&(base_gamma[5]));
      dirac_prod(&(g5_times_gQ_L[icontr]), &(base_gamma[5]),&(base_gamma[list_opQ[icontr]]));
      dirac_prod(&(gsource_times_g5_R[icontr]), &(base_gamma[gsource]),&(base_gamma[5]));
      dirac_prod(&(g5_times_gQ_R[icontr]), &(base_gamma[5]),&(base_gamma[list_opQ[icontr]]));

     //Remind that D- rotates as 1+ig5, but D-^-1 rotates as 1-ig5,
      //moreover (D^-1)^dagger rotates again as 1+ig5 (pweee!!!)

      //f1 < 1: do rotation for a quark propagator
      //f1 = 1: do nothing (physical basis already)

	// Note that in the connected part the Gamma Structure to multiply in order to rotate to physical basis is different that in the disconnected sector (ciclic trace...)
 
      if(fdL<1)
        switch(rdL)
          {
          case 0: //This is (D-^-1)^dag
            dirac_prod(&(gsource_times_g5_L[icontr]), &(gsource_times_g5_L[icontr]),&Pplus);
            dirac_prod(&(g5_times_gQ_L[icontr]), &Pplus,&(g5_times_gQ_L[icontr]));
            break;
          case 1: //This is (D+^-1)^dag
            dirac_prod(&(gsource_times_g5_L[icontr]), &(gsource_times_g5_L[icontr]),&Pminus);
            dirac_prod(&(g5_times_gQ_L[icontr]), &Pminus,&(g5_times_gQ_L[icontr]));
            break;
          }

      if(fsL<1)
        switch(rsL)
          {
          case 0: //This is D-^-1
            dirac_prod(&(g5_times_gQ_R[icontr]), &(g5_times_gQ_R[icontr]),&Pminus);
            dirac_prod(&(gsource_times_g5_L[icontr]), &Pminus,&(gsource_times_g5_L[icontr]));
            break;
          case 1: //This is D+^-1
            dirac_prod(&(g5_times_gQ_R[icontr]), &(g5_times_gQ_R[icontr]),&Pplus);
            dirac_prod(&(gsource_times_g5_L[icontr]), &Pplus,&(gsource_times_g5_L[icontr]));
            break;
          }

      if(fdR<1)
        switch(rdR)
          {
          case 0: //This is (D-^-1)^dag
            dirac_prod(&(gsource_times_g5_R[icontr]), &(gsource_times_g5_R[icontr]),&Pplus);
            dirac_prod(&(g5_times_gQ_R[icontr]), &Pplus,&(g5_times_gQ_R[icontr]));
            break;
          case 1: //This is (D+^-1)^dag
            dirac_prod(&(gsource_times_g5_R[icontr]), &(gsource_times_g5_R[icontr]),&Pminus);
            dirac_prod(&(g5_times_gQ_R[icontr]), &Pminus,&(g5_times_gQ_R[icontr]));
            break;
          }
  
      if(fsR<1)
        switch(rsR)
          {
          case 0: //This is D-^-1
            dirac_prod(&(g5_times_gQ_L[icontr]), &(g5_times_gQ_L[icontr]),&Pminus);
            dirac_prod(&(gsource_times_g5_R[icontr]), &Pminus,&(gsource_times_g5_R[icontr]));
            break;
          case 1: //This is D+^-1
            dirac_prod(&(g5_times_gQ_L[icontr]), &(g5_times_gQ_L[icontr]),&Pplus);
            dirac_prod(&(gsource_times_g5_R[icontr]), &Pplus,&(gsource_times_g5_R[icontr]));
            break;
          }

    }

 //Call for the routine which does the real contraction for the Mezzotto and the otto:
  trace_g_sdag_g_s_g_sdag_g_s(otto,gsource_times_g5_L,sdL,g5_times_gQ_L,ssL,gsource_times_g5_R,sdR,g5_times_gQ_R,ssR,ncontr);

}

int main(int narg,char **arg)
{
  int tot_prop_read=0;
  int tot_contr_made=0;
  
  double tot_reading_time=0;
  double tot_contract_time=0;

 double tic,tic1,tac,tac1;

  //Basic mpi initialization
  init_nissa();

  if(narg<2 && rank==0)
      {
	fprintf(stderr,"Use: %s input_file\n",arg[0]);
	fflush(stderr);
	MPI_Abort(MPI_COMM_WORLD,1);
      }

  open_input(arg[1]);

  ///////////////////////////////////////////////////////
  ////////////////// Basic input file ///////////////////
  ///////////////////////////////////////////////////////

  //Read the volume
  int L,T;
  read_str_int("L",&L);
  read_str_int("T",&T);

  //Read the time location of the source
  int twall;
  read_str_int("TWall",&twall);

  //Read the number of propagators of the first list
  int nprop_list1;
  read_str_int("NPropFirstlist",&nprop_list1);
  master_printf("Nprop of the first list: %d\n",nprop_list1);

  //Read the name, mass, theta and other flags for the first list
  char **base_filename1=(char**)malloc(sizeof(char*)*nprop_list1);
  for(int iprop1=0;iprop1<nprop_list1;iprop1++) base_filename1[iprop1]=(char*)malloc(1024);
  char **end_filename1=(char**)malloc(sizeof(char*)*nprop_list1);
  for(int iprop1=0;iprop1<nprop_list1;iprop1++) end_filename1[iprop1]=(char*)malloc(1024);
  double  *mass_prop1=(double*)malloc(sizeof(double)*nprop_list1);
  double *theta_prop1=(double*)malloc(sizeof(double)*nprop_list1);
  int     *phys_prop1=   (int*)malloc(sizeof(int)   *nprop_list1);
  int        *r_prop1=   (int*)malloc(sizeof(int)   *nprop_list1);
  for(int iprop=0;iprop<nprop_list1;iprop++)
    {
      read_str(base_filename1[iprop],1024);
      read_str(end_filename1[iprop],1024);
      read_double(&(mass_prop1[iprop]));
      read_double(&(theta_prop1[iprop]));
      read_int(&(phys_prop1[iprop]));
      read_int(&(r_prop1[iprop]));

      if(debug_lvl && rank==0)
	printf(" prop.%d %s %s, m=%f th=%f phys=%d r=%d\n",iprop,base_filename1[iprop],end_filename1[iprop],mass_prop1[iprop],theta_prop1[iprop],phys_prop1[iprop],r_prop1[iprop]);
    }
      
  //Read the number of propagators of the second list
  int nprop_list2;
  read_str_int("NPropSecondlist",&nprop_list2);
  if(rank==0) printf("Nprop of the second list: %d\n",nprop_list2);

  //Read the name, mass, theta and other flags for the second list
  char **base_filename2=(char**)malloc(sizeof(char*)*nprop_list2);
  for(int iprop2=0;iprop2<nprop_list2;iprop2++) base_filename2[iprop2]=(char*)malloc(1024);
  char **end_filename2=(char**)malloc(sizeof(char*)*nprop_list2);
  for(int iprop2=0;iprop2<nprop_list2;iprop2++) end_filename2[iprop2]=(char*)malloc(1024);
  double  *mass_prop2=(double*)malloc(sizeof(double)*nprop_list2);
  double *theta_prop2=(double*)malloc(sizeof(double)*nprop_list2);
  int     *phys_prop2=   (int*)malloc(sizeof(int)   *nprop_list2);
  int        *r_prop2=   (int*)malloc(sizeof(int)   *nprop_list2);
  for(int iprop=0;iprop<nprop_list2;iprop++)
    {
      read_str(base_filename2[iprop],1024);
      read_str(end_filename2[iprop],1024);
      read_double(&(mass_prop2[iprop]));
      read_double(&(theta_prop2[iprop]));
      read_int(&(phys_prop2[iprop]));
      read_int(&(r_prop2[iprop]));

      if(debug_lvl && rank==0)
	printf(" prop.%d %s %s, m=%f th=%f phys=%d r=%d\n",iprop,base_filename2[iprop],end_filename2[iprop],mass_prop2[iprop],theta_prop2[iprop],phys_prop2[iprop],r_prop2[iprop]);
    }

  //Read the number of propagators of the third list
  int nprop_list3;
  read_str_int("NPropThirdlist",&nprop_list3);
  master_printf("Nprop of the third list: %d\n",nprop_list3);

  //Read the name, mass, theta and other flags for the third list
  char **base_filename3=(char**)malloc(sizeof(char*)*nprop_list3);
  for(int iprop3=0;iprop3<nprop_list3;iprop3++) base_filename3[iprop3]=(char*)malloc(1024);
  char **end_filename3=(char**)malloc(sizeof(char*)*nprop_list3);
  for(int iprop3=0;iprop3<nprop_list3;iprop3++) end_filename3[iprop3]=(char*)malloc(1024);
  double  *mass_prop3=(double*)malloc(sizeof(double)*nprop_list3);
  double *theta_prop3=(double*)malloc(sizeof(double)*nprop_list3);
  int     *phys_prop3=   (int*)malloc(sizeof(int)   *nprop_list3);
  int        *r_prop3=   (int*)malloc(sizeof(int)   *nprop_list3);
  for(int iprop=0;iprop<nprop_list3;iprop++)
    {
      read_str(base_filename3[iprop],1024);
      read_str(end_filename3[iprop],1024);
      read_double(&(mass_prop3[iprop]));
      read_double(&(theta_prop3[iprop]));
      read_int(&(phys_prop3[iprop]));
      read_int(&(r_prop3[iprop]));

      if(debug_lvl && rank==0)
        printf(" prop.%d %s %s, m=%f th=%f phys=%d r=%d\n",iprop,base_filename3[iprop],end_filename3[iprop],mass_prop3[iprop],theta_prop3[iprop],phys_prop3[iprop],r_prop3[iprop]);
    }

  //Read the number of propagators of the fourth list
  int nprop_list4;
  read_str_int("NPropFourthlist",&nprop_list4);
  if(rank==0) printf("Nprop of the fourth list: %d\n",nprop_list4);

  //Read the name, mass, theta and other flags for the fourth list
  char **base_filename4=(char**)malloc(sizeof(char*)*nprop_list4);
  for(int iprop4=0;iprop4<nprop_list4;iprop4++) base_filename4[iprop4]=(char*)malloc(1024);
  char **end_filename4=(char**)malloc(sizeof(char*)*nprop_list4);
  for(int iprop4=0;iprop4<nprop_list4;iprop4++) end_filename4[iprop4]=(char*)malloc(1024);
  double  *mass_prop4=(double*)malloc(sizeof(double)*nprop_list4);
  double *theta_prop4=(double*)malloc(sizeof(double)*nprop_list4);
  int     *phys_prop4=   (int*)malloc(sizeof(int)   *nprop_list4);
  int        *r_prop4=   (int*)malloc(sizeof(int)   *nprop_list4);
  for(int iprop=0;iprop<nprop_list4;iprop++)
    {
      read_str(base_filename4[iprop],1024);
      read_str(end_filename4[iprop],1024);
      read_double(&(mass_prop4[iprop]));
      read_double(&(theta_prop4[iprop]));
      read_int(&(phys_prop4[iprop]));
      read_int(&(r_prop4[iprop]));

      if(debug_lvl && rank==0)
        printf(" prop.%d %s %s, m=%f th=%f phys=%d r=%d\n",iprop,base_filename4[iprop],end_filename4[iprop],mass_prop4[iprop],theta_prop4[iprop],phys_prop4[iprop],r_prop4[iprop]);
    }

  //Check the number of propagators in the lists
  if(nprop_list1 != nprop_list3)
    {
      if(rank==0)
	{
	  fprintf(stderr,"Error: the first and third lists must have the same number of propagators\n");
	  fflush(stderr);
	  MPI_Abort(MPI_COMM_WORLD,1);
	}
    }

  if(nprop_list2 != nprop_list4)
    {
      if(rank==0)
	{
	  fprintf(stderr,"Error: the second and foruth lists must have the same number of propagators\n");
	  fflush(stderr);
	  MPI_Abort(MPI_COMM_WORLD,1);
	}
    }


  //Read the number of contractions

  int ncontr;
  read_str_int("NContr_eight",&ncontr);
  if(rank==0) printf("Number of contractions for the eight: %d\n",ncontr);

  //Initialize the list of correlations and the list of operators
  //contiguous allocation for the three points
  complex **mezzottos=(complex**)malloc(sizeof(complex*)*ncontr);
  mezzottos[0]=(complex*)malloc(sizeof(complex)*ncontr*glb_size[0]); 
  complex **ottos=(complex**)malloc(sizeof(complex*)*ncontr);
  ottos[0]=(complex*)malloc(sizeof(complex)*ncontr*glb_size[0]);

  int *op=(int*)malloc(sizeof(int)*ncontr);

  for(int icontr=0;icontr<ncontr;icontr++)
    {
      mezzottos[icontr]=mezzottos[0]+icontr*glb_size[0];
      ottos[icontr]=ottos[0]+icontr*glb_size[0];

      read_int(&(op[icontr]));

      if(rank==0 && debug_lvl) printf(" eight contr.%d %d \n",icontr,op[icontr]);
    }

  int ncontr2;
  read_str_int("NContr_2points",&ncontr2);
  if(rank==0) printf("Number of contractions for the 2points: %d\n",ncontr2);


  complex **mezzottoR=(complex**)malloc(sizeof(complex*)*ncontr2);
  mezzottoR[0]=(complex*)malloc(sizeof(complex)*ncontr2*glb_size[0]);

  complex **mezzottoL=(complex**)malloc(sizeof(complex*)*ncontr2);
  mezzottoL[0]=(complex*)malloc(sizeof(complex)*ncontr2*glb_size[0]);

  int *op1=(int*)malloc(sizeof(int)*ncontr2);
  int *op2=(int*)malloc(sizeof(int)*ncontr2);
  for(int icontr=0;icontr<ncontr2;icontr++)
    {

      mezzottoL[icontr]=mezzottoL[0]+icontr*glb_size[0];
      mezzottoR[icontr]=mezzottoR[0]+icontr*glb_size[0];

      //Read the operator pairs
      read_int(&(op1[icontr]));
      read_int(&(op2[icontr]));

      if(rank==0 && debug_lvl) printf(" 2points contr.%d %d %d\n",icontr,op1[icontr],op2[icontr]);
    }

  //Read the output filename for the tree-points
  char outfile[1024];
  read_str_str("Output_eight",outfile,1024);

  //Read the output filename for the two-points
  char outfile2[1024];
  read_str_str("Output_2points",outfile2,1024);

  close_input();

  //Init the MPI grid 
  init_grid(T,L);
  
  //Calculate the number of blocks for the first list
  int nprop_per_block=compute_allocable_propagators_list(nprop_list1);
  int nblocks=nprop_list1/(nprop_per_block);
  if(nprop_list1>nblocks*nprop_per_block) nblocks++;

  //allocate the spinors
  colorspinspin **spinor1=(colorspinspin**)malloc(sizeof(colorspinspin*)*nprop_per_block);
  for(int iprop1=0;iprop1<nprop_per_block;iprop1++) spinor1[iprop1]=(colorspinspin*)malloc(sizeof(colorspinspin)*loc_vol);
  colorspinspin **spinor3=(colorspinspin**)malloc(sizeof(colorspinspin*)*nprop_per_block);
  for(int iprop3=0;iprop3<nprop_per_block;iprop3++) spinor3[iprop3]=(colorspinspin*)malloc(sizeof(colorspinspin)*loc_vol);
  colorspinspin *spinor2=(colorspinspin*)malloc(sizeof(colorspinspin)*loc_vol);
  colorspinspin **spinor4=(colorspinspin**)malloc(sizeof(colorspinspin*)*2);
  for(int r4=0;r4<2;r4++) spinor4[r4]=(colorspinspin*)malloc(sizeof(colorspinspin)*loc_vol);

  colorspinspin *spinor2_ptr;    //This will point to spinor2
  colorspinspin *spinor4_ptr[2]; //This will point to spinor4

  ///////////////////////////////////////////
  
  //take initial time
  if(debug_lvl)
    {
      MPI_Barrier(cart_comm);
      tic=MPI_Wtime();
    }

  FILE *fout=NULL;
  if(rank==0)
    {
      fout=fopen(outfile,"w");
      if(fout==NULL)
	{
	  fprintf(stderr,"Couldn't open the file: %s",outfile);
	  fflush(stderr);
	  MPI_Abort(MPI_COMM_WORLD,1);
	}
    }

  FILE *fout2=NULL;
  if(rank==0)
    {
      fout2=fopen(outfile2,"w");
      if(fout2==NULL)
        {
          fprintf(stderr,"Couldn't open the file: %s",outfile2);
          fflush(stderr);
          MPI_Abort(MPI_COMM_WORLD,1);
        }
    }

  int spat_vol=glb_size[1]*glb_size[1]*glb_size[1];


// Start of the contractions

  //Loop over the blocks of the first and third lists
  for(int iblock=0;iblock<nblocks;iblock++)
    {
      int iblock_first=iblock*nprop_per_block;
      int iblock_last=min_int((iblock+1)*nprop_per_block,nprop_list1);
      int iblock_length=iblock_last-iblock_first;

      if(rank==0 && debug_lvl) printf("Block %d/%d length: %d\n",iblock+1,nblocks,iblock_length);

      //now read the whole block of the first list
      for(int iprop1=0;iprop1<iblock_length;iprop1++)
      {
        int counter1=iblock_first+iprop1;

        if(debug_lvl)
          {
            if(rank==0 && debug_lvl>1) printf("Going to read propagator %d/%d from the first list: %s\n",iprop1+1,iblock_length,base_filename1[counter1]);
            MPI_Barrier(cart_comm);
            tic1=MPI_Wtime();
          }
        read_colorspinspin(spinor1[iprop1],base_filename1[counter1],end_filename1[counter1]);
        if(debug_lvl)
          {
            MPI_Barrier(cart_comm);
            tac1=MPI_Wtime();
            tot_reading_time+=tac1-tic1;
            tot_prop_read++;
          }
      }

      //now read the whole block of the third list
      for(int iprop3=0;iprop3<iblock_length;iprop3++)
      {
        int counter3=iblock_first+iprop3;

        if(debug_lvl)
          {
            if(rank==0 && debug_lvl>1) printf("Going to read propagator %d/%d from the third list: %s\n",iprop3+1,iblock_length,base_filename1[counter3]);
            MPI_Barrier(cart_comm);
            tic1=MPI_Wtime();
          }
        read_colorspinspin(spinor3[iprop3],base_filename3[counter3],end_filename1[counter3]);
        if(debug_lvl)
          {
            MPI_Barrier(cart_comm);
            tac1=MPI_Wtime();
            tot_reading_time+=tac1-tic1;
            tot_prop_read++;
          }
      }

     //initialize variable iprop4_prev
     int iprop4_prev = -1;

     //now loop over the second propagator
     for(int iprop2=0;iprop2<nprop_list2;iprop2++)
     {
          //read the second propagator one by one
          spinor2_ptr=spinor2;
 	  // look for this propagator in the first list (LEFT SIDE)
          for(int iprop1=0;iprop1<iblock_length;iprop1++)
            {
              int counter1=iblock_first+iprop1;
              if(strcmp(base_filename1[counter1],base_filename2[iprop2])==0 && strcmp(end_filename1[counter1],end_filename2[iprop2])==0)
                {
                  spinor2_ptr=spinor1[iprop1];
                  if(debug_lvl && rank==0) printf("Propagator %s found in the position %d of the first list\n",base_filename2[iprop2],counter1);
                }
            }
          //if not found in the first list, load it
          if(spinor2_ptr==spinor2)
            {
              if(debug_lvl)
                {
                  if(rank==0) printf("Going to read propagator %d/%d from the second list: %s %s \n",iprop2+1,nprop_list2,base_filename2[iprop2], end_filename2[iprop2]);
                  MPI_Barrier(cart_comm);
                  tic1=MPI_Wtime();
                }
              read_colorspinspin(spinor2,base_filename2[iprop2],end_filename2[iprop2]);
              if(debug_lvl)
                {
                  MPI_Barrier(cart_comm);
                  tac1=MPI_Wtime();
                  tot_reading_time+=tac1-tic1;
                  tot_prop_read++;
                }
            }

           int iprop4 = iprop2 - r_prop2[iprop2];
          //read the two flavors for the fourth propagator
           if(iprop4 != iprop4_prev)
           {
                  spinor4_ptr[0]=spinor4[0];
                  spinor4_ptr[1]=spinor4[1];
		  // look for this propagator among the ones of the third list (RIGHT SIDE)
          for(int iprop3=0;iprop3<iblock_length;iprop3+=2)
            {
              int counter3=iblock_first+iprop3;
              if(strcmp(base_filename3[counter3],base_filename4[iprop4])==0 && strcmp(end_filename3[counter3],end_filename4[iprop4])==0)
                {
                  spinor4_ptr[0]=spinor3[iprop3];
                  if(debug_lvl && rank==0) printf("Propagator %s found in the position %d of the third list\n",base_filename4[iprop4],counter3);
                  spinor4_ptr[1]=spinor3[iprop3+1];
                  if(debug_lvl && rank==0) printf("Propagator %s found in the position %d of the third list\n",base_filename4[iprop4+1],counter3+1);
                }
            }
              //if not found among the third ones, load it
	          if(spinor4_ptr[0]==spinor4[0])
        	    {
                  for(int r4=0;r4<2;r4++)
                  {
                     int counter4=iprop4+r4;

	             	 if(debug_lvl)
        	        {
               		  if(rank==0) printf("Going to read propagator %d/%d from the fourth list: %s %s \n",counter4+1,nprop_list4,base_filename4[counter4], end_filename4[counter4]);
	                  MPI_Barrier(cart_comm);
        		  tic1=MPI_Wtime();
               		 }

	              read_colorspinspin(spinor4[r4],base_filename4[counter4],end_filename4[counter4]);
	              if(debug_lvl)
        	        {
                	  MPI_Barrier(cart_comm);
	                  tac1=MPI_Wtime();
        	          tot_reading_time+=tac1-tic1;
                	  tot_prop_read++;
               		 }
               	  }
           	    }
           	 iprop4_prev=iprop4;
           }


		      for(int iprop1=0;iprop1<iblock_length;iprop1++)
			{
			  int counter1=iblock_first+iprop1;

		//Calculate two points (LEFT SIDE)
		      if (rank==0) printf ("Computing 2 points LEFT... \n");

			  if(rank==0)
			    fprintf(fout2," # m1=%f th1=%f r1=%d , m2=%f th2=%f r2=%d\n",mass_prop1[counter1],theta_prop1[counter1],r_prop1[counter1],mass_prop2[iprop2],theta_prop2[iprop2],r_prop2[iprop2]);
			  if(debug_lvl)
			    {
			      if(rank==0 && debug_lvl>1) printf("Going to perform (prop%d,prop%d) contractions\n",iprop1+1,iprop2+1);
			      MPI_Barrier(cart_comm);
			      tic1=MPI_Wtime();
			    }
             		 meson_two_points(mezzottoL,op1,spinor1[iprop1],op2,spinor2_ptr,ncontr2,phys_prop1[counter1],r_prop1[counter1],phys_prop2[iprop2],r_prop2[iprop2]);
		         if(debug_lvl)
               		 {
			   MPI_Barrier(cart_comm);
			   tac1=MPI_Wtime();
			   tot_contract_time+=tac1-tic1;
	                  tot_contr_made+=ncontr2;
			 }
			 
			 if(rank==0)
			   {
			     
			     for(int icontr=0;icontr<ncontr2;icontr++)
			       {
				 fprintf(fout2,"\n");
				 fprintf(fout2," # %s%s\n",gtag[op1[icontr]],gtag[op2[icontr]]);
				 for(int tempt=0;tempt<glb_size[0];tempt++)
        	         {
                	          int t=tempt+twall;
                        	  if(t>=glb_size[0]) t-=glb_size[0];

				  fprintf(fout2,"%+016.16g\t%+016.16g\n",mezzottoL[icontr][t][0]/spat_vol,mezzottoL[icontr][t][1]/spat_vol);
				     }
			      }
                  fprintf(fout2,"\n");
                  if(debug_lvl>1) fflush(fout2);
			   }

    //fix the index of the third propagator with flavor r3=0 and with the same mass of the current one from the first list
          int iprop3_base = iprop1 - r_prop1[counter1];

     //now loop over the two flavors of the third propagator
          for(int r3=0;r3<2;r3++)
          {
              int iprop3=iprop3_base+r3;
              int counter3=iblock_first+iprop3;

	 	  //now select the fourth propagator with the same mass as the current one from the second list and with the appropriate flavor
                int r123 = r_prop1[counter1] + r_prop2[iprop2] + r_prop3[counter3];
                int r4 =  1 - r123 + 2 * (r123 / 2);
                int counter4 = iprop4 + r4;


		//Calculate two points (RIGHT SIDE)
               if (rank==0) printf ("Computing 2 points RIGHT... \n");

                         if(rank==0)
                         fprintf(fout2," # m3=%f th3=%f r3=%d , m4=%f th4=%f r4=%d\n",mass_prop3[counter3],theta_prop3[counter3],r_prop3[counter3],mass_prop4[counter4],theta_prop4[counter4],r_prop4[counter4]);
                         if(debug_lvl)
                         {
                          if(rank==0 && debug_lvl>1) printf("Going to perform (prop%d,prop%d) contractions\n",counter3+1,counter4+1);
                          MPI_Barrier(cart_comm);
                          tic1=MPI_Wtime();
                         }
                         meson_two_points(mezzottoR,op1,spinor3[iprop3],op2,spinor4_ptr[r4],ncontr2,phys_prop3[counter3],r_prop3[counter3],phys_prop4[counter4],r_prop4[counter4]);
                         if(debug_lvl)
                         {
                          MPI_Barrier(cart_comm);
                          tac1=MPI_Wtime();
                          tot_contract_time+=tac1-tic1;
                          tot_contr_made+=ncontr2;
                         }

                        if(rank==0)
                        {
                          for(int icontr=0;icontr<ncontr2;icontr++)
                          {
                              fprintf(fout2,"\n");
                              fprintf(fout2," # %s%s\n",gtag[op1[icontr]],gtag[op2[icontr]]);
                              for(int tempt=0;tempt<glb_size[0];tempt++)
                                {
                                  int t=tempt+twall;
                                  if(t>=glb_size[0]) t-=glb_size[0];

                                 fprintf(fout2,"%+016.16g\t%+016.16g\n",mezzottoR[icontr][t][0]/spat_vol,mezzottoR[icontr][t][1]/spat_vol);
                                 }
                          }
                             fprintf(fout2,"\n");
                             if(debug_lvl>1) fflush(fout2);
                      }

                 //Calculate all the three points 

		                   if(rank==0)
		                   {
		                   printf ("Computing 3 points ... \n");
                		   fprintf(fout," # m1=%f r1=%d , m2=%f r2=%d , m3=%f r3=%d , m4=%f r4=%d\n",mass_prop1[counter1],r_prop1[counter1],mass_prop2[iprop2],r_prop2[iprop2],mass_prop3[counter3],r_prop3[counter3],mass_prop4[counter4],r_prop4[counter4]);
	                       }
                		   if(debug_lvl)
                  		   {
                    			   if(rank==0) printf("Going to perform (prop%d,prop%d,prop%d,prop%d) contractions\n",counter1+1,iprop2+1,counter3+1,counter4+1);
                			   MPI_Barrier(cart_comm);
                      			   tic1=MPI_Wtime();
                  		   }

                		  BK_eight_disconnected(mezzottos,op,spinor1[iprop1],spinor2_ptr,spinor3[iprop3],spinor4_ptr[r4],ncontr,phys_prop1[counter1],r_prop1[counter1],phys_prop2[iprop2],r_prop2[iprop2],phys_prop3[counter3],r_prop3[counter3],phys_prop4[counter4],r_prop4[counter4]);
                   		  BK_eight_connected(ottos,op,spinor1[iprop1],spinor2_ptr,spinor3[iprop3],spinor4_ptr[r4],ncontr,phys_prop1[counter1],r_prop1[counter1],phys_prop2[iprop2],r_prop2[iprop2],phys_prop3[counter3],r_prop3[counter3],phys_prop4[counter4],r_prop4[counter4]);
		                  if(debug_lvl)
                  		 {
                		     MPI_Barrier(cart_comm);
                 		     tac1=MPI_Wtime();
		                     tot_contract_time+=tac1-tic1;
                		     tot_contr_made+=2*ncontr;
                  		 }

               			if(rank==0)
                  		{
                    		 for(int icontr=0;icontr<ncontr;icontr++)
                      		 {
                        	   fprintf(fout,"\n");
                  		       fprintf(fout," # %s%s-%s%s DISCONNECTED \n",gtag[5],gtag[op[icontr]],gtag[5],gtag[op[icontr]]);
		                       for(int tempt=0;tempt<glb_size[0];tempt++)
                         		{
		                            int t=tempt+twall;
                		            if(t>=glb_size[0]) t-=glb_size[0];

                        		    fprintf(fout,"%+016.16g\t%+016.16g\n",mezzottos[icontr][t][0]/spat_vol,mezzottos[icontr][t][1]/spat_vol);
                         		 }
                        	       fprintf(fout,"\n");
                        	       fprintf(fout," # %s%s-%s%s CONNECTED\n",gtag[5],gtag[op[icontr]],gtag[5],gtag[op[icontr]]);
                       		       for(int tempt=0;tempt<glb_size[0];tempt++)
                    		       {
                       			     int t=tempt+twall;
                        		     if(t>=glb_size[0]) t-=glb_size[0];

		                             fprintf(fout,"%+016.16g\t%+016.16g\n",ottos[icontr][t][0]/spat_vol,ottos[icontr][t][1]/spat_vol);
                         		  }

				               }
                		   fprintf(fout,"\n");
                 		   if(debug_lvl>1) fflush(fout);
		                }

     } // end for over r3
    } // end for over iprop1

   } // end for over iprop2

  } // end for over iblock

  //take final time
  if(debug_lvl)
    {
      MPI_Barrier(cart_comm);
      tac=MPI_Wtime();
      if(rank==0)
	{
	  printf("\nTotal time elapsed: %f s of which:\n",tac-tic);
	  printf(" - %f s (%2.2f/100) to read %d propagators  (aver. %f s/prop) s\n",
		 tot_reading_time,tot_reading_time/(tac-tic)*100,tot_prop_read,tot_reading_time/tot_prop_read);
	  printf(" - %f s (%2.2f/100) to make %d contractions (aver. %f s/contr) s\n",
		 tot_contract_time,tot_contract_time/(tac-tic)*100,tot_contr_made,tot_contract_time/tot_contr_made);
	}
    }

  ///////////////////////////////////////////

  if(rank==0)
    {
      printf("\nEverything ok, exiting!\n");
      fclose(fout);
      fclose(fout2);
    }

  /*
  free(mass_prop4);
  free(theta_prop4);
  free(phys_prop4);
  free(r_prop4);
  for(int r4=0;r4<2;r4++) free(spinor4[r4]);
  free(spinor4);
  for(int iprop4=0;iprop4<nprop_list4;iprop4++) free(base_filename4[iprop4]);
  free(base_filename4);

  free(mass_prop2);
  free(theta_prop2);
  free(phys_prop2);
  free(r_prop2);
  free(spinor2);
  for(int iprop2=0;iprop2<nprop_list2;iprop2++) free(base_filename2[iprop2]);
  free(base_filename2);

  free(mass_prop3);
  free(theta_prop3);
  free(phys_prop3);
  free(r_prop3);
  for(int iprop3=0;iprop3<nprop_per_block;iprop3++) free(spinor3[iprop3]);
  free(spinor3);
  for(int iprop3=0;iprop3<nprop_list3;iprop3++) free(base_filename3[iprop3]);
  free(base_filename3);

  free(mass_prop1);
  free(theta_prop1);
  free(phys_prop1);
  free(r_prop1);
  for(int iprop1=0;iprop1<nprop_per_block;iprop1++) free(spinor1[iprop1]);
  free(spinor1);
  for(int iprop1=0;iprop1<nprop_list1;iprop1++) free(base_filename1[iprop1]);
  free(base_filename1);

  free(mezzottos[0]);
  free(mezzottos);

  free(ottos[0]);
  free(ottos);

  free(mezzottoL[0]);
  free(mezzottoL);

  free(mezzottoR[0]);
  free(mezzottoR);

  free(op);
  */

  close_nissa();

  return 0;
}
