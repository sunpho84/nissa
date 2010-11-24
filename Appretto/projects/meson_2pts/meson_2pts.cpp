
//This program calculate a list of two points function contracting all
//the propagators in the first list with all the propagators in the
//second list. The propagators on the second list are loaded one by
//one. The first list is split into block, each of them as large as
//much it is possible.

#include <mpi.h>
#include <fstream>
#include <lemon.h>

#include "appretto.h"

using namespace std;

//Calculate the number of blocks in which to divide the first list
int compute_nblocks_first_list(int nprop_list1)
{
  //First of all it checks if there is enough room for two propagators.
  //This is the minimal requirement for the program to be able to work.
  colorspinspin *fuf;
  fuf=(colorspinspin*)malloc(2*sizeof(colorspinspin)*loc_vol);

  if(fuf==NULL)
    {
      if(rank==0) cerr<<"Error: not enough memory for two propagators"<<endl;
      MPI_Abort(MPI_COMM_WORLD,1);
    }
  else if(debug and rank==0) cout<<"Ok there is enough memory to load 2 propagators"<<endl;

  free(fuf);

  //Now determine the largest number of propagator of the first list (and one additional) loadable at once.
  //We are sure that we can allocate at least 2 props, so it will exit with at least output=1.
  int nblocks=0;
  int nprop_per_block;
  do
    {
      nblocks++;
      nprop_per_block=nprop_list1/nblocks;
      fuf=(colorspinspin*)malloc((nprop_per_block+1)*sizeof(colorspinspin)*loc_vol);
    }
  while(fuf==NULL);

  if(debug and rank==0)
    cout<<"Will split the "<<nprop_list1<<" propagators in the first list in "<<nblocks<<" blocks of "<<nprop_per_block<<" each"<<endl;

  return nblocks;
}

//Wrapper to the real calculation
//This function takes care to make the revert on the SECOND spinor, putting both the needed gamma5
//It also put the rotators to the physical basis if asked
void meson_two_points(complex *corr,int *list_op1,colorspinspin *s1,int *list_op2,colorspinspin *s2,int ncontr,int f1,int r1,int f2,int r2)
{
  //Temporary vectors for the internal gamma
  dirac_matr t1[ncontr],t2[ncontr];

  for(int icontr=0;icontr<ncontr;icontr++)
    {
      
      //Put the two gamma5 needed for the revert of the second spinor
      dirac_prod(t1[icontr], base_gamma[5],base_gamma[list_op1[icontr]]);
      dirac_prod(t2[icontr], base_gamma[list_op2[icontr]],base_gamma[5]);
      
      //Remember that D- rotate as 1+ig5, but D-^-1 rotate ad 1-ig5,
      //morover (D^-1)^dagger rotate again as 1+ig5 (pweee!!!)

      if(f1!=1)
	switch(r1)
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

      if(f2!=1)
	switch(r2)
	  {
	  case 0: //This is (D-^-1)^dagger
	    dirac_prod(t2[icontr], t2[icontr],Pplus);
	    dirac_prod(t1[icontr], Pplus,t1[icontr]);
	    break;
	  case 1: //This is (D+^-1)^dagger
	    dirac_prod(t2[icontr], t2[icontr],Pminus);
	    dirac_prod(t1[icontr], Pminus,t1[icontr]);
	    break;
	  }
    }
  
  //Call for the routine which do the real contraction
  trace_g_s_g_sdag(corr,t1,s1,t2,s2,ncontr);
}

int main(int narg,char **arg)
{
  //Basic mpi initialization
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

  //Read the position of the source
  int twall;
  read_int("TWall",twall);

  //Read the number of propagator of the first list
  int nprop_list1;
  read_int("NPropFirstlist",nprop_list1);
  if(rank==0) cout<<"Nprop of the first list: "<<nprop_list1<<endl;

  //Read the name, mass etc of the first list
  char base_filename1[nprop_list1][1024];
  double mass_prop1[nprop_list1];
  double theta_prop1[nprop_list1];
  int phys_prop1[nprop_list1];
  int r_prop1[nprop_list1];
  for(int iprop=0;iprop<nprop_list1;iprop++)
    {
      read_str(base_filename1[iprop],1024);
      read_double(mass_prop1[iprop]);
      read_double(theta_prop1[iprop]);
      read_int(phys_prop1[iprop]);
      read_int(r_prop1[iprop]);

      if(debug and rank==0)
	cout<<" prop "<<iprop<<", m="<<mass_prop1[iprop]<<" th="<<theta_prop1[iprop]<<" phys="<<phys_prop1[iprop]<<" r="<<r_prop1[iprop]<<endl;
    }
      
  //Read the number of propagator of the first list
  int nprop_list2;
  read_int("NPropSecondlist",nprop_list2);
  if(rank==0) cout<<"Nprop of the first list: "<<nprop_list2<<endl;

  //Read the name, mass etc of the first list
  char base_filename2[nprop_list2][1024];
  double mass_prop2[nprop_list2];
  double theta_prop2[nprop_list2];
  int phys_prop2[nprop_list2];
  int r_prop2[nprop_list2];
  for(int iprop=0;iprop<nprop_list2;iprop++)
    {
      read_str(base_filename2[iprop],1024);
      read_double(mass_prop2[iprop]);
      read_double(theta_prop2[iprop]);
      read_int(phys_prop2[iprop]);
      read_int(r_prop2[iprop]);

      if(debug and rank==0)
	cout<<" prop "<<iprop<<", m="<<mass_prop2[iprop]<<" th="<<theta_prop2[iprop]<<" phys="<<phys_prop2[iprop]<<" r="<<r_prop2[iprop]<<endl;
    }
      
  //Read the number of contractions
  int ncontr;
  read_int("NContr",ncontr);

  //Initialize the list of correlations and the list of operators
  complex contr[ncontr][glb_size[0]];
  int *op1=new int[ncontr]; 
  int *op2=new int[ncontr]; 
  for(int icontr=0;icontr<ncontr;icontr++)
    {
      //Read the operators pairs
      read_int(op1[icontr]);
      read_int(op2[icontr]);

      if(rank==0) cout<<" contr. "<<icontr<<" "<<op1[icontr]<<" "<<op2[icontr]<<endl;
    }

  close_input();

  //Init the MPI grid 
  init_grid();
  
  //Calculate the number of blocks in which to divide the first list
  int nblocks=compute_nblock_first_list(nprop_list1);
  int nprop_per_block=nprop_list1/nblocks;
  if(nprop_list1>nblock*nprop_per_block) nblocks++;

  //allocate the spinors
  colorspinspin *spinor1[nprop_per_block];
  for(int iprop1=0;iprop1<nprop_per_block;iprop1++) spinor1[iprop1]=new colorspinspin[loc_vol];
  colorspinspin *spinor2=new colorspinspin[loc_vol];

  ///////////////////////////////////////////
  
  //take initial time
  double tic;
  if(debug)
    {
      MPI_Barrier(cart_comm);
      tic=MPI_Wtime();
    }

  ofstream fout;
  fout.precision(12);
  fout.width(12);
  if(rank==0) fout.open("two_points_correlations");
  int spat_vol=glb_size[1]*glb_size[1]*glb_size[1];

  //Loop over the blocks of the first list
  for(int iblock=0;iblock<nblocks;iblock++)
    {
      int iblock_first=iblock*nprop_per_block;
      int iblock_last=min((iblock+1)*nprop_per_block,nprop_list1);
      int iblock_length=iblock_last-iblock_first;

      //now read the whole first block
      for(int iprop1=0;iprop1<iblock_length;iprop1++)
      {
	int counter=iblock_first+iprop1;
	
	read_colorspinspin(base_filename1[counter],spinor1[iprop1]);
      }

      //now loop over the second popagator
      for(int iprop2=0;iprop2<nprop_list2;iprop2++)
	{
	  //read the second propagator one by one
	  read_colorspinspin(base_filename2[iprop2],spinor2);
	  
	  //Calculate all the two points between spinor 1 and spinor2
	  for(int iprop1=0;iprop1<iblock_length;iprop1++)
	    {
	      int counter=iblock_first+iprop1;

	      if(rank==0)
		fout<<" #"
		    <<" m1="<<mass_prop1[counter]<<" th1="<<theta_prop1[counter]<<" r1="<<r_prop1[counter]<<" ,"
		    <<" m2="<<mass_prop2[iprop2]<<" th2="<<theta_prop2[iprop2]<<" r2="<<r_prop2[iprop2]<<endl<<endl;

	      meson_two_points((complex*)contr,op1,spinor1[iprop1],op2,spinor2,ncontr,phys_prop1[counter],r_prop1[counter],phys_prop2[iprop2],r_prop2[iprop2]);

	      if(rank==0)
		{
		  fout<<endl;
		  
		  for(int icontr=0;icontr<ncontr;icontr++)
		    {
		      fout<<" # "<<op1[icontr]<<" "<<op2[icontr]<<endl;
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
	}
    }

  //take final time
  double tac;
  if(debug)
    {
      MPI_Barrier(cart_comm);
      tac=MPI_Wtime();
      if(rank==0) cout<<"Time elapsed in doing "<<nprop_list1*nprop_list2*ncontr<<" contraction: "<<tac-tic<<" s"<<endl;
    }

  ///////////////////////////////////////////

  for(int iprop1=0;iprop1<nprop_per_block;iprop1++) delete[] spinor1[iprop1];
  delete[] spinor2;

  delete[] op1;
  delete[] op2;

  close_appretto();

  return 0;
}
