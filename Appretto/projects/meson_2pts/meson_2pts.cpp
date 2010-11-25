
//This program calculates a list of two-point functions contracting 
//all the propagators in the first list with all the propagators in
//the second list. The propagators of the second list are loaded one 
//by one. The first list is split into blocks, each of them as large 
//as possible.

#include <mpi.h>
#include <fstream>
#include <lemon.h>

#include "appretto.h"

using namespace std;

//Calculate the maximum number of allocable propagators
int compute_allocable_propagators(int nprop_list)
{
  //First of all check if there is enough room for two propagators.
  //This is the minimal requirement for the program to be able to work.
  colorspinspin *fuf;
  fuf=(colorspinspin*)malloc(2*sizeof(colorspinspin)*loc_vol);

  if(fuf==NULL)
    {
      if(rank==0) cerr<<"Error: not enough memory for two propagators"<<endl;
      MPI_Abort(MPI_COMM_WORLD,1);
    }
  else if(debug and rank==0) cout<<"Ok there is enough memory to load two propagators"<<endl;

  free(fuf);

  //Now determine the largest number of propagator of the first list (and one additional) loadable at once.
  //We are sure that we can allocate at least 2 props, so it will exit with at least nprop_max=1.
  int nprop_max=nprop_list + 1;
  do
    {
      nprop_max--;
      fuf=(colorspinspin*)malloc((nprop_max+1)*sizeof(colorspinspin)*loc_vol);
    }
  while(fuf==NULL);

  free(fuf);

  if(debug and rank==0)
    cout<<"Will allocate "<<nprop_max<<" propagators from a list with "<<nprop_list<<" propagators"<<endl;

  return nprop_max;
}

//This function takes care to make the revert on the SECOND spinor, putting the needed gamma5
//It also applies the appropriate rotators to the physical basis if asked
void meson_two_points(complex *corr,int *list_op1,colorspinspin *s1,int *list_op2,colorspinspin *s2,int ncontr,int f1,int r1,int f2,int r2)
{
  //Temporary vectors for the internal gamma
  dirac_matr t1[ncontr],t2[ncontr];

  for(int icontr=0;icontr<ncontr;icontr++)
    {
      
      //Put the two gamma5 needed for the revert of the second spinor
      dirac_prod(t1[icontr], base_gamma[5],base_gamma[list_op1[icontr]]);
      dirac_prod(t2[icontr], base_gamma[list_op2[icontr]],base_gamma[5]);
      
      //Remind that D- rotates as 1+ig5, but D-^-1 rotates as 1-ig5,
      //moreover (D^-1)^dagger rotates again as 1+ig5 (pweee!!!)

      //f1 < 1: do rotation for a quark propagator
      //f1 = 1: do nothing (physical basis already)
      //f1 > 1: do rotation for a (charged) sequential propagator

      if(f1<1)
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

      if(f1>1)
        switch(r1)
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

      if(f2<1)
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

     if(f2>1)
        switch(r2)
          {
          case 0: //This is (D-^-1)^dagger
            dirac_prod(t2[icontr], t2[icontr],Pplus);
            dirac_prod(t1[icontr], Pminus,t1[icontr]);
            break;
          case 1: //This is (D+^-1)^dagger
            dirac_prod(t2[icontr], t2[icontr],Pminus);
            dirac_prod(t1[icontr], Pplus,t1[icontr]);
            break;
          }

    }
  
  //Call for the routine which does the real contraction
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

  //Read the time location of the source
  int twall;
  read_int("TWall",twall);

  //Read the number of propagators of the first list
  int nprop_list1;
  read_int("NPropFirstlist",nprop_list1);
  if(rank==0) cout<<"Nprop of the first list: "<<nprop_list1<<endl;

  //Read the name, mass, theta and other flags for the first list
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
      
  //Read the number of propagators of the second list
  int nprop_list2;
  read_int("NPropSecondlist",nprop_list2);
  if(rank==0) cout<<"Nprop of the first list: "<<nprop_list2<<endl;

  //Read the name, mass, theta and other flags for the second list
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
      //Read the operator pairs
      read_int(op1[icontr]);
      read_int(op2[icontr]);

      if(rank==0) cout<<" contr. "<<icontr<<" "<<op1[icontr]<<" "<<op2[icontr]<<endl;
    }

  close_input();

  //Init the MPI grid 
  init_grid();
  
  //Calculate the number of blocks for the first list
  int nprop_per_block=compute_allocable_propagators(nprop_list1);
  int nblocks=nprop_list1/nprop_per_block;
  if(nprop_list1>nblocks*nprop_per_block) nblocks++;

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
  fout.precision(16);
  fout.width(16);
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
	  colorspinspin *spinor2_ptr; //This will point to spinor2 if the prop. is not in the first list
	  //check if the file is already loaded
	  spinor2_ptr=spinor2;
	  for(int iprop1=0;iprop1<iblock_length;iprop1++)
	    {
	      int counter=iblock_first+iprop1;
	      if(strcmp(base_filename1[counter],base_filename2[iprop2])==0)
		{
		  spinor2_ptr=spinor1[iprop1];
		  if(debug and rank==0) cout<<"Propagator "<<base_filename2[iprop2]<<" found in the position "<<counter<<" of the first list"<<endl;
		}
	    }
	  //if not found in the first list, load it
	  if(spinor2_ptr==spinor2) read_colorspinspin(base_filename2[iprop2],spinor2);
	  
	  //Calculate all the two points between spinor 1 and spinor2
	  for(int iprop1=0;iprop1<iblock_length;iprop1++)
	    {
	      int counter=iblock_first+iprop1;

	      if(rank==0)
		fout<<" #"
		    <<" m1="<<mass_prop1[counter]<<" th1="<<theta_prop1[counter]<<" r1="<<r_prop1[counter]<<" ,"
		    <<" m2="<<mass_prop2[iprop2] <<" th2="<<theta_prop2[iprop2] <<" r2="<<r_prop2[iprop2]<<endl;

	      meson_two_points((complex*)contr,op1,spinor1[iprop1],op2,spinor2_ptr,ncontr,phys_prop1[counter],r_prop1[counter],phys_prop2[iprop2],r_prop2[iprop2]);

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

 if(rank==0) fout.close();

  for(int iprop1=0;iprop1<nprop_per_block;iprop1++) delete[] spinor1[iprop1];
  delete[] spinor2;

  delete[] op1;
  delete[] op2;

  close_appretto();

  return 0;
}
