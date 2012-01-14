
//This program calculates a list of two-point functions contracting 
//all the propagators in the first list with all the propagators in
//the second list. The propagators of the second list are loaded one 
//by one. The first list is split into blocks, each of them as large 
//as possible.

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

//Calculate the maximum number of allocable propagators
//First of all check if there is enough room for the configuration,
//then for the two/three propagators.
//This is the minimal requirement for the program to be able to work.
int compute_allocable_propagators(int nprop_list,int nch_contr)
{
  quad_su3 *temp_conf=NULL;
  as2t_su3 *temp_clov=NULL;
  if(nch_contr>0)
    {
      temp_conf=(quad_su3*)malloc(sizeof(quad_su3)*(loc_vol+loc_bord+loc_edge));
      if(temp_conf==NULL && rank>0)
	{
	  fprintf(stderr,"Unable to allocate the space for the gauge configuration!\n");
	  MPI_Abort(MPI_COMM_WORLD,1);
	}

      temp_clov=(as2t_su3*)malloc(sizeof(as2t_su3)*loc_vol);
      if(temp_clov==NULL && rank>0)
	{
	  fprintf(stderr,"Unable to allocate the space for the P_munu term!\n");
	  MPI_Abort(MPI_COMM_WORLD,1);
	}
    }

  colorspinspin *fuf=NULL;
  int nmin_req;
  if(nch_contr==0) nmin_req=2;
  else nmin_req=3;
  fuf=(colorspinspin*)malloc(nmin_req*sizeof(colorspinspin)*loc_vol);

  if(fuf==NULL)
    {
      if(rank==0)
	{
	  fprintf(stderr,"Error: not enough memory for %d propagators\n",nmin_req);
	  fflush(stderr);
	  MPI_Abort(MPI_COMM_WORLD,1);
	}
    }
  else if(debug_lvl>1 && rank==0) printf("Ok there is enough memory to load %d propagators\n",nmin_req);

  free(fuf);

  //Now determine the largest number of propagator of the first list (and one additional) loadable at once.
  //We are sure that we can allocate at least nmin_req props, so it will exit with at least nprop_max=1.
  int nprop_max=nprop_list+nmin_req-1;
  do
    {
      nprop_max--;
      fuf=(colorspinspin*)malloc((nprop_max+1)*sizeof(colorspinspin)*loc_vol);
    }
  while(fuf==NULL);

  free(fuf);

  if(debug_lvl && rank==0)
    printf("Will allocate %d propagators from a list with %d propagators\n",nprop_max,nprop_list);

  if(nch_contr>0)
    {
      free(temp_conf);
      free(temp_clov);
    }

  return nprop_max;
}

//This function takes care to make the revert on the FIRST spinor, putting the needed gamma5
//It also applies the appropriate rotators to the physical basis if asked
void meson_two_points(complex *corr,int *list_op1,colorspinspin *s1,int *list_op2,colorspinspin *s2,int ncontr,int f1,int r1,int f2,int r2)
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
  trace_g_sdag_g_s(corr,t1,s1,t2,s2,ncontr);
}

int main(int narg,char **arg)
{
  int tot_prop_read=0;
  int tot_contr_made=0;
  
  double tot_reading_time=0;
  double tot_contract_time=0;

  //Basic mpi initialization
  init_nissa();

  if(narg<2 && rank==0)
      {
	fprintf(stderr,"Use: %s input_file\n",arg[0]);
	fflush(stderr);
	MPI_Abort(MPI_COMM_WORLD,1);
      }

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
  int    * phys_prop1=   (int*)malloc(sizeof(int)   *nprop_list1);
  int        *r_prop1=   (int*)malloc(sizeof(int)   *nprop_list1);
  for(int iprop=0;iprop<nprop_list1;iprop++)
    {
      read_str(base_filename1[iprop],1024);
      read_double(&(mass_prop1[iprop]));
      read_double(&(theta_prop1[iprop]));
      read_int(&(phys_prop1[iprop]));
      read_int(&(r_prop1[iprop]));

      if(debug_lvl && rank==0)
	printf(" prop.%d %s, m=%f th=%f phys=%d r=%d\n",iprop,base_filename1[iprop],mass_prop1[iprop],theta_prop1[iprop],phys_prop1[iprop],r_prop1[iprop]);
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
  int    * phys_prop2=   (int*)malloc(sizeof(int)   *nprop_list2);
  int        *r_prop2=   (int*)malloc(sizeof(int)   *nprop_list2);
  for(int iprop=0;iprop<nprop_list2;iprop++)
    {
      read_str(base_filename2[iprop],1024);
      read_double(&(mass_prop2[iprop]));
      read_double(&(theta_prop2[iprop]));
      read_int(&(phys_prop2[iprop]));
      read_int(&(r_prop2[iprop]));

      if(debug_lvl && rank==0)
	printf(" prop.%d %s, m=%f th=%f phys=%d r=%d\n",iprop,base_filename2[iprop],mass_prop2[iprop],theta_prop2[iprop],phys_prop2[iprop],r_prop2[iprop]);
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
  
  //Read the location of the gauge configuration if needed
  char gaugeconf_file[1024];
  if(nch_contr>0) read_str_str("GaugeConf",gaugeconf_file,1024);
    
  //Read the output filename
  char outfile[1024];
  read_str_str("Output",outfile,1024);
  
  close_input();
  
  //Init the MPI grid 
  init_grid(T,L);
  
  //Calculate the number of blocks for the first list
  int nprop_per_block=compute_allocable_propagators(nprop_list1,nch_contr);
  int nblocks=nprop_list1/nprop_per_block;
  if(nprop_list1>nblocks*nprop_per_block) nblocks++;
  
  //allocate the spinors
  colorspinspin **spinor1=(colorspinspin**)malloc(sizeof(colorspinspin*)*nprop_per_block);
  for(int iprop1=0;iprop1<nprop_per_block;iprop1++) spinor1[iprop1]=(colorspinspin*)malloc(sizeof(colorspinspin)*loc_vol);
  colorspinspin *spinor2=(colorspinspin*)malloc(sizeof(colorspinspin)*loc_vol);
  
  //if we have to calculate the chromo-magnetic operator allocate one additional spinor
  //if necessary allocate and load the gauge configuration,and allocate the space for the pmunu term
  colorspinspin *ch_spinor=NULL;
  quad_su3 *gauge_conf=NULL;
  as2t_su3 *Pmunu=NULL;
  if(nch_contr!=0)
    {
      ch_spinor=(colorspinspin*)malloc(sizeof(colorspinspin)*loc_vol);
      gauge_conf=(quad_su3*)malloc(sizeof(quad_su3)*(loc_vol+loc_bord+loc_edge));
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
  
  //if necessary, load the gauge configuration and calculate the pmunu term
  if(nch_contr>0)
    {
      read_ildg_gauge_conf(gauge_conf,gaugeconf_file);
      communicate_lx_gauge_borders(gauge_conf);
      communicate_lx_gauge_edges(gauge_conf);

      double gplaq=global_plaquette_lx_conf(gauge_conf);
      if(rank==0) printf("plaq: %.10g\n",gplaq);
      
      Pmunu_term(Pmunu,gauge_conf);
      free(gauge_conf);
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
	read_colorspinspin(spinor1[iprop1],base_filename1[counter],NULL);
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
	  colorspinspin *spinor2_ptr; //This will point to spinor2 if the prop. is not in the first list
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
	      read_colorspinspin(spinor2,base_filename2[iprop2],NULL);
	      if(debug_lvl)
		{
		  MPI_Barrier(cart_comm);
		  tac1=MPI_Wtime();
		  tot_reading_time+=tac1-tic1;
		  tot_prop_read++;
		}
	    }
	  
	  //apply the chromo magnetic operator to the second spinor
	  if(nch_contr>0) unsafe_apply_chromo_operator_to_colorspinspin(ch_spinor,Pmunu,spinor2_ptr);
	  
	  //Calculate all the two points between spinor 1 and spinor2
	  for(int iprop1=0;iprop1<iblock_length;iprop1++)
	    {
	      int counter=iblock_first+iprop1;

	      if(rank==0)
		fprintf(fout," # m1=%f th1=%f r1=%d , m2=%f th2=%f r2=%d\n",mass_prop1[counter],theta_prop1[counter],r_prop1[counter],mass_prop2[iprop2],theta_prop2[iprop2],r_prop2[iprop2]);

	      double tic1=0,tac1;
	      if(debug_lvl)
		{
		  if(rank==0 && debug_lvl>1) printf("Going to perform (prop%d,prop%d) contractions\n",iprop1+1,iprop2+1);
		  MPI_Barrier(cart_comm);
		  tic1=MPI_Wtime();
		}
	      meson_two_points(contr,op1,spinor1[iprop1],op2,spinor2_ptr,ncontr,phys_prop1[counter],r_prop1[counter],phys_prop2[iprop2],r_prop2[iprop2]);
	      if(nch_contr>0) meson_two_points(ch_contr,ch_op1,spinor1[iprop1],ch_op2,ch_spinor,nch_contr,phys_prop1[counter],r_prop1[counter],phys_prop2[iprop2],r_prop2[iprop2]);

	      if(debug_lvl)
		{
		  MPI_Barrier(cart_comm);
		  tac1=MPI_Wtime();
		  tot_contract_time+=tac1-tic1;
		  tot_contr_made+=ncontr+nch_contr;
		}

	      if(rank==0)
		{
		  print_contractions_to_file(fout,ncontr,op1,op2,contr,twall,"");
		  if(nch_contr>0) print_contractions_to_file(fout,nch_contr,ch_op1,ch_op2,ch_contr,twall,"CHROMO-");

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
  free(phys_prop2);
  free(r_prop2);
  free(spinor2);
  for(int iprop2=0;iprop2<nprop_list2;iprop2++) free(base_filename2[iprop2]);
  free(base_filename2);

  free(mass_prop1);
  free(theta_prop1);
  free(phys_prop1);
  free(r_prop1);
  for(int iprop1=0;iprop1<nprop_per_block;iprop1++) free(spinor1[iprop1]);
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
