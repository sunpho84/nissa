// This program create the two sources, invert them and perform all the needed contractions

//////////////////////////////////////////////////////////////////////////////////////////////////////
// Three-point disconnected part (Mezzotos)
//               _                                          _      
//   Sum_x{ Tr [ S(md;x,tL) gsource S(ms,x,tL) OPL  ] Tr [ S(md;x,tR) gsource S(ms,x,tR) OP_R  ]  }      
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
                .....                           .....
              ..     ..                       ..     ..
             .         .         (x)         .         .        
   gsource  X           X  OP_L        OP_R  X           X gsource
     (tL)    .         .               (tR)  .         .
              ..     ..                       ..      ..
                .....                           .....
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

char conf_path[1024];
quad_su3 *conf;
double kappa;

int nmass;
double *mass;
double put_theta[4],old_theta[4]={0,0,0,0};

//source data
int seed,noise_type,twall[2];
double rad2=1.414213562373095048801688724209;
spincolor *source;
colorspinspin *original_source;

//vectors for the spinor data
colorspinspin **S[2][2];
spincolor **cgmms_solution,*temp_vec[2];

//cgmms inverter parameters
double stopping_residue;
double minimal_residue;
int stopping_criterion;
int niter_max;

//two points contractions
complex **contr_otto,**contr_mezzotto,*contr_2pts;
char outfile_otto[1024],outfile_mezzotto[1024],outfile_2pts[1024];
int *op1_2pts,*op2_2pts;
int ncontr_2pts;

//timings
int ninv_tot=0,ncontr_tot=0;
double tot_time=0,inv_time=0,contr_time=0;

//number of spectator masses
int nspec;

//Generate the source for the dirac index
void generate_source(int LR)
{ //reset
  memset(original_source,0,sizeof(colorspinspin)*loc_vol);
 
  for(int loc_site=0;loc_site<loc_vol;loc_site++)
    if(glb_coord_of_loclx[loc_site][0]==twall[LR])
      for(int ic=0;ic<3;ic++)
        { //real part
          if(noise_type>=2) original_source[loc_site][ic][0][0][0]=pm_one(loc_site)/rad2;
          else original_source[loc_site][ic][0][0][0]=noise_type;
          //imaginary part
          if(noise_type==4) original_source[loc_site][ic][0][0][1]=pm_one(loc_site)/rad2;
         
          for(int d=1;d<4;d++) //copy the other three dirac indexes
            memcpy(original_source[loc_site][ic][d][d],original_source[loc_site][ic][0][0],sizeof(complex));
        }
}

//Parse all the input file
void initialize_Bk(char *input_path)
{
  open_input(input_path);

  // 1) Read information about the gauge conf
 
  //Read the volume
  int L,T;
  read_str_int("L",&L);
  read_str_int("T",&T);
  //Init the MPI grid
  init_grid(T,L);
  //Gauge path
  read_str_str("GaugeConfPath",conf_path,1024);
  //Kappa
  read_str_double("Kappa",&kappa);
 
  // 2) Read information about the source
 
  //Read the seed and initialize the random generator
  read_str_int("Seed",&seed);
  init_random(seed);
  //Read the location of the left wall
  read_str_int("TLeftWall",&(twall[0]));
  read_str_int("TRightWall",&(twall[1]));
  //Read the noise type
  read_str_int("NoiseType",&noise_type);

  // 3) Read list of masses

  read_list_of_doubles("NMass",&nmass,&mass);

  // 4) Info about the inverter

  //Residue
  read_str_double("Residue",&stopping_residue);
  //Stopping criterion
  stopping_criterion=numb_known_stopping_criterion;
  char str_stopping_criterion[1024];
  read_str_str("StoppingCriterion",str_stopping_criterion,1024);
  int isc=0;
  do
    {
      if(strcasecmp(list_known_stopping_criterion[isc],str_stopping_criterion)==0) stopping_criterion=isc;
      isc++;
    }
  while(isc<numb_known_stopping_criterion && stopping_criterion==numb_known_stopping_criterion);
 
  if(stopping_criterion==numb_known_stopping_criterion && rank==0)
    {
      fprintf(stderr,"Unknown stopping criterion: %s\n",str_stopping_criterion);
      fprintf(stderr,"List of known stopping criterions:\n");
      for(int isc=0;isc<numb_known_stopping_criterion;isc++) fprintf(stderr," %s\n",list_known_stopping_criterion[isc]);
      MPI_Abort(MPI_COMM_WORLD,1);
    }
  if(stopping_criterion==sc_standard) read_str_double("MinimalResidue",&minimal_residue);
 
  //Number of iterations
  read_str_int("NiterMax",&niter_max);
 
  // 5) contraction list for eight

  contr_otto=(complex**)malloc(sizeof(complex*)*16);
  contr_otto[0]=(complex*)malloc(sizeof(complex)*16*glb_size[0]); 
  contr_mezzotto=(complex**)malloc(sizeof(complex*)*16);
  contr_mezzotto[0]=(complex*)malloc(sizeof(complex)*16*glb_size[0]); 
 
  for(int iop=0;iop<16;iop++)
    {
      contr_otto[iop]=contr_otto[0]+iop*glb_size[0];
      contr_mezzotto[iop]=contr_mezzotto[0]+iop*glb_size[0];
    }
  
  read_str_str("OutfileOtto",outfile_otto,1024);
  read_str_str("OutfileTwoPoints",outfile_2pts,1024);
  read_str_int("NSpec",&nspec);
  read_str_int("NContrTwoPoints",&ncontr_2pts);
  contr_2pts=(complex*)malloc(sizeof(complex)*ncontr_2pts*glb_size[0]);
  op1_2pts=(int*)malloc(sizeof(int)*ncontr_2pts);
  op2_2pts=(int*)malloc(sizeof(int)*ncontr_2pts);
  for(int icontr=0;icontr<ncontr_2pts;icontr++)
    {

      //Read the operator pairs
      read_int(&(op1_2pts[icontr]));
      read_int(&(op2_2pts[icontr]));

      if(debug_lvl) master_printf(" contr.%d %d %d\n",icontr,op1_2pts[icontr],op2_2pts[icontr]);
    }


  close_input();

  ////////////////////////////////////// end of input reading/////////////////////////////////

  //allocate gauge conf, Pmunu and all the needed spincolor and colorspinspin
  conf=nissa_malloc("conf",loc_vol+loc_bord,quad_su3);

  //load the gauge conf, propagate borders, calculate plaquette and PmuNu term
  read_ildg_gauge_conf(conf,conf_path);
  communicate_lx_gauge_borders(conf);
 
  double gplaq=global_plaquette_lx_conf(conf);
  master_printf("plaq: %.18g\n",gplaq);
 
  //Put the anti-periodic condition on the temporal border
  put_theta[0]=1;
  adapt_theta(conf,old_theta,put_theta,1,0);

  //Allocate all the propagators colorspinspin vectors
  for(int LR=0;LR<2;LR++)
    for(int UD=0;UD<2;UD++)
      {
        S[LR][UD]=(colorspinspin**)malloc(sizeof(colorspinspin*)*nmass);
        for(int iprop=0;iprop<nmass;iprop++) S[LR][UD][iprop]=nissa_malloc("S",loc_vol,colorspinspin);
      }
 
  //Allocate nmass spincolors, for the cgmms solutions
  cgmms_solution=(spincolor**)malloc(sizeof(spincolor*)*nmass);
  for(int imass=0;imass<nmass;imass++) cgmms_solution[imass]=nissa_malloc("cgmms_solution",loc_vol+loc_bord,spincolor);
  temp_vec[0]=nissa_malloc("temp_vec[0]",loc_vol,spincolor);
  temp_vec[1]=nissa_malloc("temp_vec[1]",loc_vol,spincolor);
 
  //Allocate one spincolor for the source
  source=nissa_malloc("source",loc_vol+loc_bord,spincolor);

  //Allocate the colorspinspin of the LR source
  original_source=nissa_malloc("original_source",loc_vol,colorspinspin);
}

//Finalization
void close_Bk()
{
  if(rank==0)
    {
      printf("\n");
      printf("Total time: %g, of which:\n",tot_time);
      printf(" - %02.2f%s to perform %d inversions (%2.2gs avg)\n",inv_time/tot_time*100,"%",ninv_tot,inv_time/ninv_tot);
      printf(" - %02.2f%s to perform %d contr. (%2.2gs avg)\n",contr_time/tot_time*100,"%",ncontr_tot,contr_time/ncontr_tot);
    }

  nissa_free(conf);
  for(int iprop=0;iprop<nmass;iprop++) for(int LR=0;LR<2;LR++) for(int UD=0;UD<2;UD++) nissa_free(S[LR][UD][iprop]);
  nissa_free(temp_vec[0]);nissa_free(temp_vec[1]);

  close_nissa();
}

//calculate the propagators
void calculate_S(int LR)
{
  for(int id=0;id<4;id++)
    { //loop over the source dirac index
      for(int ivol=0;ivol<loc_vol;ivol++)
        {
          get_spincolor_from_colorspinspin(source[ivol],original_source[ivol],id);
          //put the g5
          for(int id1=2;id1<4;id1++) for(int ic=0;ic<3;ic++) for(int ri=0;ri<2;ri++) source[ivol][id1][ic][ri]*=-1;
        }
     
      communicate_lx_spincolor_borders(source);
     
      double part_time=-take_time();
      inv_Q2_cgmms(cgmms_solution,source,NULL,conf,kappa,mass,nmass,niter_max,stopping_residue,minimal_residue,stopping_criterion);
      part_time+=take_time();ninv_tot++;inv_time+=part_time;
      char tag[2][10]={"left","right"};
      master_printf("Finished the %s source inversion, dirac index %d in %g sec\n",tag[LR],id,part_time);
     
      for(int imass=0;imass<nmass;imass++)
        { //reconstruct the doublet
          reconstruct_doublet(temp_vec[0],temp_vec[1],cgmms_solution[imass],conf,kappa,mass[imass]);
          master_printf("Mass %d (%g) reconstructed \n",imass,mass[imass]);
          for(int r=0;r<2;r++) //convert the id-th spincolor into the colorspinspin
            for(int i=0;i<loc_vol;i++)
              put_spincolor_into_colorspinspin(S[LR][r][imass][i],temp_vec[r][i],id);
        }
    }
 
  for(int r=0;r<2;r++) //remember that D^-1 rotate opposite than D!
    for(int ipropS=0;ipropS<nmass;ipropS++) //put the (1+ig5)/sqrt(2) factor
      rotate_vol_colorspinspin_to_physical_basis(S[LR][r][ipropS],!r,!r);
}

void Bk_eights(colorspinspin *SL1,colorspinspin *SL2,colorspinspin *SR1,colorspinspin *SR2)
{
  //Temporary vectors for the internal gamma
  dirac_matr tsource[16],tsink[16];

  for(int igamma=0;igamma<16;igamma++)
    {
      //Put the two gamma5 needed for the revert of the d spinor
      tsource[igamma]=base_gamma[0]; //g5*g5
      dirac_prod(&(tsink[igamma]),&(base_gamma[5]),&(base_gamma[igamma]));
    }
 
  //Call the routine which does the real contraction for the Mezzotto and the Otto
  trace_g_sdag_g_s_g_sdag_g_s(contr_otto,tsource,SL1,tsink,SL2,tsource,SR1,tsink,SR2,16);
  sum_trace_g_sdag_g_s_times_trace_g_sdag_g_s(contr_mezzotto,tsource,SL1,tsink,SL2,tsource,SR1,tsink,SR2,16);
}

void meson_two_points(colorspinspin *s1,colorspinspin *s2)
{
  //Temporary vectors for the internal gamma
  dirac_matr t1[ncontr_2pts],t2[ncontr_2pts];
 
  for(int icontr=0;icontr<ncontr_2pts;icontr++)
    {
      //Put the two gamma5 needed for the revert of the first spinor
      dirac_prod(&(t1[icontr]),&(base_gamma[op1_2pts[icontr]]),&(base_gamma[5]));
      dirac_prod(&(t2[icontr]),&(base_gamma[5]),&(base_gamma[op2_2pts[icontr]]));
    }
  //Call for the routine which does the real contraction
  trace_g_sdag_g_s(contr_2pts,t1,s1,t2,s2,ncontr_2pts);
}

//print all the passed contractions to the file
void print_ottos_contractions_to_file(FILE *fout)
{
  double norm=glb_size[1]*glb_size[2]*glb_size[3];
 
  if(rank==0)
    for(int icontr=0;icontr<16;icontr++)
      {
        fprintf(fout,"\n");
        print_contraction_to_file(fout,icontr,5,contr_mezzotto[icontr],twall[0],"DISCONNECTED ",norm);
        fprintf(fout,"\n");
        print_contraction_to_file(fout,icontr,5,contr_otto[icontr],twall[0],"CONNECTED ",norm);
      }
}

//print all the passed contractions to the file
void print_two_points_contractions_to_file(FILE *fout,int LR)
{
  double norm=glb_size[1]*glb_size[2]*glb_size[3];
 
  if(rank==0)
    for(int icontr=0;icontr<ncontr_2pts;icontr++)
      {
        fprintf(fout,"\n");
        print_contraction_to_file(fout,op2_2pts[icontr],op1_2pts[icontr],contr_2pts+icontr*glb_size[0],twall[0],"",norm);
      }
}

//Calculate and print to file all the contractions
void calculate_all_contractions()
{
  FILE *fout_otto=open_text_file_for_output(outfile_otto);
  FILE *fout_2pts=open_text_file_for_output(outfile_2pts);
 
  contr_time-=take_time();

  //Ottos
  for(int im2=0;im2<nmass;im2++)
   for(int r2=0;r2<2;r2++)
    for(int im1=0;im1<nspec;im1++)
     for(int r1=0;r1<2;r1++)
          for(int r3=0;r3<2;r3++)
            {
              int r4=1-(r1+r2+r3)%2;
           
              if(rank==0) fprintf(fout_otto," # m1=%f r1=%d , m2=%f r2=%d , m3=%f r3=%d , m4=%f r4=%d\n",
                                  mass[im1],r1,mass[im2],r2,mass[im1],r3,mass[im2],r4);

              Bk_eights(S[0][r1][im1],S[0][r2][im2],S[1][r3][im1],S[1][r4][im2]);
              ncontr_tot+=32;
              if(rank==0)
                {
                  print_ottos_contractions_to_file(fout_otto);
                  fprintf(fout_otto,"\n");
                }
            }

  //two points
  char tag_LR[2][10]={"LEFT","RIGHT"};
  for(int im2=0;im2<nmass;im2++)
   for(int r2=0;r2<2;r2++)
    for(int im1=0;im1<nspec;im1++)
     for(int r1=0;r1<2;r1++)
      for(int LR=0;LR<2;LR++)
            {
              if(rank==0) fprintf(fout_2pts," # m1=%f r1=%d , m2=%f r2=%d %s\n",
                                  mass[im1],r1,mass[im2],r2,tag_LR[LR]);
              meson_two_points(S[LR][r1][im1],S[LR][r2][im2]);
              if(rank==0)
                {
                  print_two_points_contractions_to_file(fout_2pts,LR);
                  fprintf(fout_2pts,"\n");
                }
              ncontr_tot+=ncontr_2pts;
            }

  contr_time+=take_time();
 
  if(rank==0)
    {
      fclose(fout_otto);
      fclose(fout_2pts);
    }
}

int main(int narg,char **arg)
{
  //Basic mpi initialization
  init_nissa();

  if(narg<2 && rank==0)
      {
        fprintf(stderr,"Use: %s input_file\n",arg[0]);
        fflush(stderr);
        MPI_Abort(MPI_COMM_WORLD,1);
      }

  tot_time-=take_time();
  initialize_Bk(arg[1]);
 
  for(int LR=0;LR<2;LR++)
    {
      generate_source(LR);
      calculate_S(LR);
    }

  calculate_all_contractions();
 
  tot_time+=take_time();
  close_Bk();
 
  return 0;
}


