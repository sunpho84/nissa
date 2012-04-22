
//This program calculates a list of two-point functions contracting 
//all the propagators in the first list with all the propagators in
//the second list. The propagators of the second list are loaded one 
//by one. The first list is split into blocks, each of them as large 
//as possible.

//After being loaded, all propagators are reconstructed, that is,
//splitted

//First list has to be the shortest.
//The meson has quark content:
//            _
//            q(m1)q(m2)
//
//This is the reference scheme:

/*               +              
                S(m1)                           
                .....          
              ..     ..        
             .         .       
        op1  X           X  op2
             .         .       
              ..     ..        
                .....          
                S(m2)          
                                    
source |------>---->----->---->| sink

*/



#include <mpi.h>
#include <lemon.h>

#include "nissa.h"
#include "common.cpp"

int main(int narg,char **arg)
{
  int tot_prop_read=0;
  int tot_contr_made=0;
  
  double gauge_reading_time=0;
  double tot_reading_time=0;
  double tot_contract_time=0;

  //services variable for reconstruction
  double theta[4]={1,0,0,0};
  double old_theta[4]={0,0,0,0};
  
  //Basic mpi initialization
  init_nissa();

  if(narg<2) crash("Use: %s input_file",arg[0]);
  
  open_input(arg[1]);

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
  if(rank==0) printf("Nprop of the first list: %d\n",nprop_list1);

  //Read the name, mass, theta and other flags for the first list
  char **base_filename1=(char**)malloc(sizeof(char*)*nprop_list1);
  for(int iprop1=0;iprop1<nprop_list1;iprop1++) base_filename1[iprop1]=(char*)malloc(1024);
  double  *mass_prop1=(double*)malloc(sizeof(double)*nprop_list1);
  double *theta_prop1=(double*)malloc(sizeof(double)*nprop_list1);
  for(int iprop=0;iprop<nprop_list1;iprop++)
    {
      read_str(base_filename1[iprop],1024);
      read_double(&(mass_prop1[iprop]));
      read_double(&(theta_prop1[iprop]));

      if(debug_lvl && rank==0)
	printf(" prop.%d %s, m=%f th=%f\n",iprop,base_filename1[iprop],mass_prop1[iprop],theta_prop1[iprop]);
    }
      
  //Read the number of propagators of the second list
  int nprop_list2;
  read_str_int("NPropSecondlist",&nprop_list2);
  if(rank==0) printf("Nprop of the second list: %d\n",nprop_list2);

  //Read the name, mass, theta and other flags for the second list
  char **base_filename2=(char**)malloc(sizeof(char*)*nprop_list2);
  for(int iprop2=0;iprop2<nprop_list2;iprop2++) base_filename2[iprop2]=(char*)malloc(1024);
  double  *mass_prop2=(double*)malloc(sizeof(double)*nprop_list2);
  double *theta_prop2=(double*)malloc(sizeof(double)*nprop_list2);
  for(int iprop=0;iprop<nprop_list2;iprop++)
    {
      read_str(base_filename2[iprop],1024);
      read_double(&(mass_prop2[iprop]));
      read_double(&(theta_prop2[iprop]));

      if(debug_lvl && rank==0)
	printf(" prop.%d %s, m=%f th=%f\n",iprop,base_filename2[iprop],mass_prop2[iprop],theta_prop2[iprop]);
    }
      
  //Read the number of contractions
  int ncontr;
  read_str_int("NContr",&ncontr);
  if(rank==0) printf("Number of contractions: %d\n",ncontr);

  //Initialize the list of correlations and the list of operators
  //contiguous allocation
  complex *contr=(complex*)malloc(sizeof(complex)*ncontr*glb_size[0]); 
  int *op1=(int*)malloc(sizeof(int)*ncontr);
  int *op2=(int*)malloc(sizeof(int)*ncontr);
  for(int icontr=0;icontr<ncontr;icontr++)
    {
      //Read the operator pairs
      read_int(&(op1[icontr]));
      read_int(&(op2[icontr]));

      if(rank==0 && debug_lvl) printf(" contr.%d %d %d\n",icontr,op1[icontr],op2[icontr]);
    }
  
  //Read the number of contractions with insertion of the chromo-magnetic operator
  int nch_contr;
  read_str_int("NChromoContr",&nch_contr);
  if(rank==0) printf("Number of chromo-contractions: %d\n",nch_contr);
  
  //Initialize the list of chromo correlations and the list of operators
  //contiguous allocation
  complex *ch_contr=(complex*)malloc(sizeof(complex)*nch_contr*glb_size[0]); 
  int *ch_op1=(int*)malloc(sizeof(int)*nch_contr);
  int *ch_op2=(int*)malloc(sizeof(int)*nch_contr);
  for(int ich_contr=0;ich_contr<nch_contr;ich_contr++)
    {
      //Read the operator pairs
      read_int(&(ch_op1[ich_contr]));
      read_int(&(ch_op2[ich_contr]));
      
      if(rank==0 && debug_lvl) printf(" chromo contr.%d %d %d\n",ich_contr,ch_op1[ich_contr],ch_op2[ich_contr]);
    }
  
  //Read the location of the gauge configuration
  char gaugeconf_file[1024];
  read_str_str("GaugeConf",gaugeconf_file,1024);

  //Read kappa
  double kappa;
  read_str_double("Kappa",&kappa);
    
  //Read the output filename
  char outfile[1024];
  read_str_str("Output",outfile,1024);
  
  close_input();
  
  //Init the MPI grid 
  init_grid(T,L);
  
  //Calculate the number of blocks for the first list
  int nprop_per_block=compute_allocable_propagators(nprop_list1,nch_contr,5);
  int nblocks=nprop_list1/nprop_per_block;
  if(nprop_list1>nblocks*nprop_per_block) nblocks++;
  
  //allocate the spinors
  colorspinspin ***spinor1=(colorspinspin***)malloc(sizeof(colorspinspin**)*nprop_per_block);
  for(int iprop1=0;iprop1<nprop_per_block;iprop1++)
    {
      spinor1[iprop1]=(colorspinspin**)malloc(sizeof(colorspinspin*)*2);
      spinor1[iprop1][0]=(colorspinspin*)malloc(sizeof(colorspinspin)*loc_vol);
      spinor1[iprop1][1]=(colorspinspin*)malloc(sizeof(colorspinspin)*loc_vol);
    }
  colorspinspin **spinor2=(colorspinspin**)malloc(sizeof(colorspinspin*)*2);
  spinor2[0]=(colorspinspin*)malloc(sizeof(colorspinspin)*loc_vol);
  spinor2[1]=(colorspinspin*)malloc(sizeof(colorspinspin)*loc_vol);

  //if we have to calculate the chromo-magnetic operator allocate one additional spinor
  //if necessary allocate and load the gauge configuration,and allocate the space for the pmunu term
  colorspinspin **ch_spinor=NULL;
  quad_su3 *gauge_conf=(quad_su3*)malloc(sizeof(quad_su3)*(loc_vol+bord_vol+edge_vol));
  as2t_su3 *Pmunu=NULL;
  if(nch_contr!=0)
    {
      ch_spinor=(colorspinspin**)malloc(sizeof(colorspinspin*)*2);
      ch_spinor[0]=(colorspinspin*)malloc(sizeof(colorspinspin)*loc_vol);
      ch_spinor[1]=(colorspinspin*)malloc(sizeof(colorspinspin)*loc_vol);

      Pmunu=(as2t_su3*)malloc(sizeof(as2t_su3)*loc_vol);
    }  
  
  ///////////////////////////////////////////
  
  //take initial time
  double tic;
  if(debug_lvl)
    {
      MPI_Barrier(cart_comm);
      tic=MPI_Wtime();
    }
  
  //load the gauge configuration 
  read_ildg_gauge_conf(gauge_conf,gaugeconf_file);
  adapt_theta(gauge_conf,old_theta,theta,1,1);
  
  double gplaq=global_plaquette_lx_conf(gauge_conf);
  if(rank==0) printf("plaq: %.10g\n",gplaq);
  
  //if necessary calculate the pmunu term
  if(nch_contr>0) Pmunu_term(Pmunu,gauge_conf);

  if(debug_lvl)
    {
      MPI_Barrier(cart_comm);
      gauge_reading_time=MPI_Wtime()-tic;
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
  
  //Loop over the blocks of the first list
  for(int iblock=0;iblock<nblocks;iblock++)
    {
      int iblock_first=iblock*nprop_per_block;
      int iblock_last=min_int((iblock+1)*nprop_per_block,nprop_list1);
      int iblock_length=iblock_last-iblock_first;

      if(rank==0 && debug_lvl) printf("Block %d/%d length: %d\n",iblock+1,nblocks,iblock_length);

      //now read the whole first block
      for(int iprop1=0;iprop1<iblock_length;iprop1++)
      {
	int counter=iblock_first+iprop1;
	
	double tic1=0,tac1;
	if(debug_lvl)
	  {
	    if(rank==0 && debug_lvl>1) printf("Going to read propagator %d/%d: %s\n",iprop1+1,iblock_length,base_filename1[counter]);
	    MPI_Barrier(cart_comm);
	    tic1=MPI_Wtime();
	  }
	theta[1]=theta[2]=theta[3]=theta_prop1[counter];
	adapt_theta(gauge_conf,old_theta,theta,1,0);
	read_tm_colorspinspin_reconstructing(spinor1[iprop1],base_filename1[counter],NULL,gauge_conf,kappa,mass_prop2[counter]);
	if(debug_lvl)
	  {
	    MPI_Barrier(cart_comm);
	    tac1=MPI_Wtime();
	    tot_reading_time+=tac1-tic1;
	    tot_prop_read++;
	  }
      }

      //now loop over the second popagator
      for(int iprop2=0;iprop2<nprop_list2;iprop2++)
	{
	  //read the second propagator one by one
	  colorspinspin **spinor2_ptr; //This will point to spinor2 if the prop. is not in the first list
	  //check if the file is already loaded
	  spinor2_ptr=spinor2;
	  for(int iprop1=0;iprop1<iblock_length;iprop1++)
	    {
	      int counter=iblock_first+iprop1;
	      if(strcmp(base_filename1[counter],base_filename2[iprop2])==0)
		{
		  spinor2_ptr=spinor1[iprop1];
		  if(debug_lvl && rank==0) printf("Propagator %s found in the position %d of the first list\n",base_filename2[iprop2],counter);
		}
	    }
	  //if not found in the first list, load it
	  if(spinor2_ptr==spinor2)
	    {
	      double tic1=0,tac1;
	      if(debug_lvl)
		{
		  if(rank==0) printf("Going to read propagator %d/%d: %s\n",iprop2+1,nprop_list2,base_filename2[iprop2]);
		  MPI_Barrier(cart_comm);
		  tic1=MPI_Wtime();
		}
	      theta[1]=theta[2]=theta[3]=theta_prop2[iprop2];
	      adapt_theta(gauge_conf,old_theta,theta,1,0);
	      read_tm_colorspinspin_reconstructing(spinor2,base_filename2[iprop2],NULL,gauge_conf,kappa,mass_prop2[iprop2]);
	      if(debug_lvl)
		{
		  MPI_Barrier(cart_comm);
		  tac1=MPI_Wtime();
		  tot_reading_time+=tac1-tic1;
		  tot_prop_read++;
		}
	    }
	  
	  //apply the chromo magnetic operator to the second spinor
	  if(nch_contr>0)
	    {
	      unsafe_apply_chromo_operator_to_colorspinspin(ch_spinor[0],Pmunu,spinor2_ptr[0]);
	      unsafe_apply_chromo_operator_to_colorspinspin(ch_spinor[1],Pmunu,spinor2_ptr[1]);
	    }

	  //Calculate all the two points between spinor 1 and spinor2
	  for(int iprop1=0;iprop1<iblock_length;iprop1++)
	    for(int r2=0;r2<2;r2++)
	      for(int r1=0;r1<2;r1++)
		{
		  int counter=iblock_first+iprop1;
		  
		  if(rank==0)
		    fprintf(fout," # m1=%f th1=%f r1=%d , m2=%f th2=%f r2=%d\n",mass_prop1[counter],theta_prop1[counter],r1,mass_prop2[iprop2],theta_prop2[iprop2],r2);
		  
		  double tic1=0,tac1;
		  if(debug_lvl)
		    {
		      if(rank==0 && debug_lvl>1) printf("Going to perform (prop%d,prop%d) contractions\n",iprop1+1,iprop2+1);
		      MPI_Barrier(cart_comm);
		      tic1=MPI_Wtime();
		    }
		  meson_two_points(contr,op1,spinor1[iprop1][r1],op2,spinor2_ptr[r2],ncontr,0,r1,0,r2);
		  if(nch_contr>0) meson_two_points(ch_contr,ch_op1,spinor1[iprop1][r1],ch_op2,ch_spinor[r2],nch_contr,0,r1,0,r2);
		  
		  if(debug_lvl)
		    {
		      MPI_Barrier(cart_comm);
		      tac1=MPI_Wtime();
		      tot_contract_time+=tac1-tic1;
		      tot_contr_made+=ncontr+nch_contr;
		    }
		  
		  if(rank==0)
		    {
		      print_contractions_to_file(fout,ncontr,op1,op2,contr,twall,"",1);
		      if(nch_contr>0) print_contractions_to_file(fout,nch_contr,ch_op1,ch_op2,ch_contr,twall,"CHROMO-",1);
		      
		      fprintf(fout,"\n");
		      if(debug_lvl>1) fflush(fout);
		    }
		}
	}
    }
  
  //take final time
  double tac;
  if(debug_lvl)
    {
      MPI_Barrier(cart_comm);
      tac=MPI_Wtime();
      if(rank==0)
	{
	  printf("\nTotal time elapsed: %f s of which:\n",tac-tic);
	  printf(" - %f s (%2.2f/100) to read gauge conf s\n",
		 gauge_reading_time,gauge_reading_time/(tac-tic)*100);
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
    }

  if(nch_contr!=0)
    {
      free(ch_spinor);
      free(Pmunu);
    }

  free(mass_prop2);
  free(theta_prop2);
  free(spinor2[0]);
  free(spinor2[1]);
  free(spinor2);
  for(int iprop2=0;iprop2<nprop_list2;iprop2++) free(base_filename2[iprop2]);
  free(base_filename2);

  free(mass_prop1);
  free(theta_prop1);
  for(int iprop1=0;iprop1<nprop_per_block;iprop1++)
    {
      free(spinor1[iprop1][0]);
      free(spinor1[iprop1][1]);
      free(spinor1[iprop1]);
    }
  free(spinor1);
  for(int iprop1=0;iprop1<nprop_list1;iprop1++) free(base_filename1[iprop1]);
  free(base_filename1);

  free(ch_contr);

  free(ch_op1);
  free(ch_op2);

  free(contr);

  free(op1);
  free(op2);

  close_nissa();

  return 0;
}
