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
#include <string.h>
#include <stdlib.h>

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
      temp_conf=(quad_su3*)malloc(sizeof(quad_su3)*(loc_vol+bord_vol+edge_vol));
      if(temp_conf==NULL) CRASH("Unable to allocate the space for the gauge configuration!");
      
      temp_clov=(as2t_su3*)malloc(sizeof(as2t_su3)*loc_vol);
      if(temp_clov==NULL) CRASH("Unable to allocate the space for the P_munu term!");
    }
  
  colorspinspin *fuf=NULL;
  int nmin_req;
  if(nch_contr==0) nmin_req=2;
  else nmin_req=3;
  fuf=(colorspinspin*)malloc(nmin_req*sizeof(colorspinspin)*loc_vol);
  
  if(fuf==NULL) CRASH("Error: not enough memory for %d propagators",nmin_req);
  else MASTER_PRINTF("Ok there is enough memory to load %d propagators\n",nmin_req);
  
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
  
  MASTER_PRINTF("Will allocate %d propagators from a list with %d propagators\n",nprop_max,nprop_list);
  
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
  int tot_contr=0;
  
  double tot_reading_time=0;
  double tot_contract_time=0;

  //Basic mpi initialization
  init_nissa();
  
  //Init timinig
  double tot_time=-take_time();
  
  //Open input file
  if(narg<2) CRASH("Use: %s input_file",arg[0]);
  open_input(arg[1]);
  
  //Init the MPI grid 
  int L,T;
  read_str_int("L",&L);
  read_str_int("T",&T);
  init_grid(T,L);
  
  //Read the time location of the source
  int twall;
  read_str_int("TWall",&twall);
  
  //Read the number of propagators of the first list
  int nprop_list1;
  read_str_int("NPropFirstlist",&nprop_list1);
  MASTER_PRINTF("Nprop of the first list: %d\n",nprop_list1);
  
  //Read the name, mass, theta and other flags for the first list
  char **base_filename1=nissa_malloc("base_filename1*",nprop_list1,char*);
  for(int iprop1=0;iprop1<nprop_list1;iprop1++) base_filename1[iprop1]=nissa_malloc("base_filename1",1024,char);
  double  *mass_prop1=nissa_malloc("mass_prop1",nprop_list1,double);
  double *theta_prop1=nissa_malloc("theta_prop1",nprop_list1,double);
  int    * phys_prop1=nissa_malloc("phys_prop1",nprop_list1,int);
  int        *r_prop1=nissa_malloc("r_prop1",nprop_list1,int);
  for(int iprop=0;iprop<nprop_list1;iprop++)
    {
      read_str(base_filename1[iprop],1024);
      read_double(&(mass_prop1[iprop]));
      read_double(&(theta_prop1[iprop]));
      read_int(&(phys_prop1[iprop]));
      read_int(&(r_prop1[iprop]));

      MASTER_PRINTF(" prop.%d %s, m=%f th=%f phys=%d r=%d\n",iprop,base_filename1[iprop],mass_prop1[iprop],theta_prop1[iprop],phys_prop1[iprop],r_prop1[iprop]);
    }
      
  //Read the number of propagators of the second list
  int nprop_list2;
  read_str_int("NPropSecondlist",&nprop_list2);
  MASTER_PRINTF("Nprop of the second list: %d\n",nprop_list2);
  
  //Read the name, mass, theta and other flags for the second list
  char **base_filename2=nissa_malloc("base_filename2*",nprop_list2,char*);
  for(int iprop2=0;iprop2<nprop_list2;iprop2++) base_filename2[iprop2]=nissa_malloc("base_filename2",1024,char);
  double  *mass_prop2=nissa_malloc("mass_prop2",nprop_list2,double);
  double *theta_prop2=nissa_malloc("theta_prop2",nprop_list2,double);
  int    * phys_prop2=nissa_malloc("phys_prop2",nprop_list2,int);
  int        *r_prop2=nissa_malloc("r_prop2",nprop_list2,int);
  for(int iprop=0;iprop<nprop_list2;iprop++)
    {
      read_str(base_filename2[iprop],1024);
      read_double(&(mass_prop2[iprop]));
      read_double(&(theta_prop2[iprop]));
      read_int(&(phys_prop2[iprop]));
      read_int(&(r_prop2[iprop]));

      MASTER_PRINTF(" prop.%d %s, m=%f th=%f phys=%d r=%d\n",iprop,base_filename2[iprop],mass_prop2[iprop],theta_prop2[iprop],phys_prop2[iprop],r_prop2[iprop]);
    }
      
  //Read the number of contractions
  int ncontr;
  read_str_int("NContr",&ncontr);
  MASTER_PRINTF("Number of contractions: %d\n",ncontr);

  //Initialize the list of correlations and the list of operators
  //contiguous allocation
  MASTER_PRINTF("%d\n",glb_size[0]);
  complex *contr=nissa_malloc("contr",ncontr*glb_size[0],complex);
  int *op1=nissa_malloc("op1",ncontr,int);
  int *op2=nissa_malloc("op2",ncontr,int);
  for(int icontr=0;icontr<ncontr;icontr++)
    {
      //Read the operator pairs
      read_int(&(op1[icontr]));
      read_int(&(op2[icontr]));

      MASTER_PRINTF(" contr.%d %d %d\n",icontr,op1[icontr],op2[icontr]);
    }
  
  //Read the number of contractions with insertion of the chromo-magnetic operator
  int nch_contr;
  read_str_int("NChromoContr",&nch_contr);
  MASTER_PRINTF("Number of chromo-contractions: %d\n",nch_contr);
  
  //Initialize the list of chromo correlations and the list of operators
  //contiguous allocation
  complex *ch_contr=nissa_malloc("ch_contr",nch_contr*glb_size[0],complex);
  int *ch_op1=nissa_malloc("ch_op1",nch_contr,int);
  int *ch_op2=nissa_malloc("ch_op2",nch_contr,int);
  for(int ich_contr=0;ich_contr<nch_contr;ich_contr++)
    {
      //Read the operator pairs
      read_int(&(ch_op1[ich_contr]));
      read_int(&(ch_op2[ich_contr]));
      
      MASTER_PRINTF(" chromo contr.%d %d %d\n",ich_contr,ch_op1[ich_contr],ch_op2[ich_contr]);
    }
  
  //Read the location of the gauge configuration if needed
  char gaugeconf_file[1024];
  if(nch_contr>0) read_str_str("GaugeConf",gaugeconf_file,1024);
    
  //Read the output filename
  char outfile[1024];
  read_str_str("Output",outfile,1024);
  
  close_input();
  
  /////////////////////////////////////////////////////
  
  //Calculate the number of blocks for the first list
  int nprop_per_block=compute_allocable_propagators(nprop_list1,nch_contr);
  int nblocks=nprop_list1/nprop_per_block;
  if(nprop_list1>nblocks*nprop_per_block) nblocks++;
  
  //allocate the spinors
  colorspinspin **spinor1=nissa_malloc("spinor1*",nprop_per_block,colorspinspin*);
  for(int iprop1=0;iprop1<nprop_per_block;iprop1++) spinor1[iprop1]=nissa_malloc("spinor1",loc_vol,colorspinspin);
  colorspinspin *spinor2=nissa_malloc("spinor2",loc_vol,colorspinspin);
  
  //if we have to calculate the chromo-magnetic operator allocate one additional spinor
  //if necessary allocate and load the gauge configuration,and allocate the space for the pmunu term
  colorspinspin *ch_spinor=NULL;
  as2t_su3 *Pmunu=NULL;
  quad_su3 *gauge_conf;
  if(nch_contr!=0)
    {
      ch_spinor=nissa_malloc("ch_spinor",loc_vol,colorspinspin);
      Pmunu=nissa_malloc("Pmunu",loc_vol,as2t_su3);
    }  
  
  ///////////////////////////////////////////
  
  //if necessary, load the gauge configuration and calculate the pmunu term
  if(nch_contr>0)
    {
      gauge_conf=nissa_malloc("conf",loc_vol+bord_vol+edge_vol,quad_su3);
  
      read_ildg_gauge_conf(gauge_conf,gaugeconf_file);
      
      MASTER_PRINTF("plaq: %.10g\n",global_plaquette_lx_conf(gauge_conf));
      
      Pmunu_term(Pmunu,gauge_conf);
      nissa_free(gauge_conf);
    }
  
  FILE *fout=open_text_file_for_output(outfile);
  
  //Loop over the blocks of the first list
  for(int iblock=0;iblock<nblocks;iblock++)
    {
      int iblock_first=iblock*nprop_per_block;
      int iblock_last=min_int((iblock+1)*nprop_per_block,nprop_list1);
      int iblock_length=iblock_last-iblock_first;
      
      MASTER_PRINTF("Block %d/%d length: %d\n",iblock+1,nblocks,iblock_length);
      
      //now read the whole first block
      for(int iprop1=0;iprop1<iblock_length;iprop1++)
      {
	int counter=iblock_first+iprop1;
	
	MASTER_PRINTF("Going to read propagator %d/%d: %s\n",iprop1+1,iblock_length,base_filename1[counter]);
	tot_reading_time-=take_time();
	read_colorspinspin(spinor1[iprop1],base_filename1[counter],NULL);
	tot_reading_time+=take_time();
	tot_prop_read++;
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
		  MASTER_PRINTF("Propagator %s found in the position %d of the first list\n",base_filename2[iprop2],counter);
		}
	    }
	  //if not found in the first list, load it
	  if(spinor2_ptr==spinor2)
	    {
	      MASTER_PRINTF("Going to read propagator %d/%d: %s\n",iprop2+1,nprop_list2,base_filename2[iprop2]);
	      tot_reading_time-=take_time();
	      read_colorspinspin(spinor2,base_filename2[iprop2],NULL);
	      tot_reading_time+=take_time();
	      tot_prop_read++;
	    }
	  
	  //apply the chromo magnetic operator to the second spinor
	  if(nch_contr>0) unsafe_apply_chromo_operator_to_colorspinspin(ch_spinor,Pmunu,spinor2_ptr);
	  
	  //Calculate all the two points between spinor 1 and spinor2
	  for(int iprop1=0;iprop1<iblock_length;iprop1++)
	    {
	      int counter=iblock_first+iprop1;

	      if(rank==0)
		fprintf(fout," # m1=%f th1=%f r1=%d , m2=%f th2=%f r2=%d\n",mass_prop1[counter],theta_prop1[counter],r_prop1[counter],mass_prop2[iprop2],theta_prop2[iprop2],r_prop2[iprop2]);

	      MASTER_PRINTF("Going to perform (prop%d,prop%d) contractions\n",iprop1+1,iprop2+1);
	      tot_contract_time-=take_time();
	      meson_two_points(contr,op1,spinor1[iprop1],op2,spinor2_ptr,ncontr,phys_prop1[counter],r_prop1[counter],phys_prop2[iprop2],r_prop2[iprop2]);

	      if(nch_contr>0) meson_two_points(ch_contr,ch_op1,spinor1[iprop1],ch_op2,ch_spinor,nch_contr,phys_prop1[counter],r_prop1[counter],phys_prop2[iprop2],r_prop2[iprop2]);

	      tot_contract_time+=take_time();
	      tot_contr+=ncontr+nch_contr;
	      
	      if(rank==0)
		{
		  print_contractions_to_file(fout,ncontr,op1,op2,contr,twall,"",1);
		  if(nch_contr>0) print_contractions_to_file(fout,nch_contr,ch_op1,ch_op2,ch_contr,twall,"CHROMO-",1);
		  
		  fprintf(fout,"\n");
		  if(nissa_verbosity>=3) fflush(fout);
		}
	    }
	}
    }
  
  //take final time
  tot_time+=take_time();
  MASTER_PRINTF("\nTotal time elapsed: %f s of which:\n",tot_time);
  MASTER_PRINTF(" - %f s (%2.2f/100) to read %d propagators  (aver. %f s/prop) s\n",
		tot_reading_time,tot_reading_time/tot_time*100,tot_prop_read,tot_reading_time/tot_prop_read);
  MASTER_PRINTF(" - %f s (%2.2f/100) to make %d contractions (aver. %f s/contr) s\n",
		tot_contract_time,tot_contract_time/tot_time*100,tot_contr,tot_contract_time/tot_contr);
  
  ///////////////////////////////////////////
  
  if(nch_contr!=0)
    {
      nissa_free(ch_spinor);
      nissa_free(Pmunu);
    }

  nissa_free(mass_prop2);
  nissa_free(theta_prop2);
  nissa_free(phys_prop2);
  nissa_free(r_prop2);
  nissa_free(spinor2);
  for(int iprop2=0;iprop2<nprop_list2;iprop2++) nissa_free(base_filename2[iprop2]);
  nissa_free(base_filename2);
  
  nissa_free(mass_prop1);
  nissa_free(theta_prop1);
  nissa_free(phys_prop1);
  nissa_free(r_prop1);
  for(int iprop1=0;iprop1<nprop_per_block;iprop1++) nissa_free(spinor1[iprop1]);
  nissa_free(spinor1);
  for(int iprop1=0;iprop1<nprop_list1;iprop1++) nissa_free(base_filename1[iprop1]);
  nissa_free(base_filename1);
  
  nissa_free(ch_contr);
  
  nissa_free(ch_op1);
  nissa_free(ch_op2);
  
  nissa_free(contr);
  
  nissa_free(op1);
  nissa_free(op2);
  
  MASTER_PRINTF("\nEverything ok, exiting!\n");
  
  close_nissa();
  
  return 0;
}
