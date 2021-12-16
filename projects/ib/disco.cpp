#include "nissa.hpp"

#include <random>
#include <tuple>

using namespace nissa;

/// Use or not the new generator
int useNewGenerator{true};
FieldRngStream field_rng_stream;

int nhits;
int seed;
double wall_time;

double mass;
double kappa;
double cSW;
double residue;

int useSme;
double apeAlpha,kappaSme;
int nApe,nSme;

int nanalyzed_conf;
int ngauge_conf;

quad_su3 *glb_conf;
quad_su3 *ape_conf;
CUDA_MANAGED spincolor *temp,**source_ptr,**prop_ptr;
clover_term_t *Cl;
inv_clover_term_t *invCl;

double prop_time,P5P5_time,S0P5_time,disco_time,P5P5_with_ins_time,prepare_time;

enum RunMode{INITIALIZING,SEARCHING_CONF,ANALYZING_CONF,STOP};
RunMode runMode{INITIALIZING};
const char stop_path[]="stop";

CUDA_HOST_AND_DEVICE
spincolor* &prop(const int t,const int r)
{
  return prop_ptr[r+2*t];
}

CUDA_HOST_AND_DEVICE
spincolor* &source(const int t)
{
  return source_ptr[t];
}

CUDA_MANAGED complex* loc_contr;
complex* temp_contr;

int iconf=0;
char conf_path[1024];
char outfolder[1024];
char run_file[1024];
lock_file_t<uint64_t> lock_file;
double init_time;


void init_simulation(int narg,char **arg)
{
  init_time=take_time();
  
  lock_file.init();
  
  open_input(arg[1]);
  
  int L,T;
  read_str_int("L",&L);
  read_str_int("T",&T);
  
  init_grid(T,L);
  
  read_str_double("WallTime",&wall_time);
  
  read_str_int("UseNewGenerator",&useNewGenerator);
  
  read_str_int("Seed",&seed);
  
  read_str_int("NHits",&nhits);
  
  read_str_double("Kappa",&kappa);
  
  read_str_double("cSW",&cSW);
  
  read_str_double("Mass",&mass);
  
  read_str_double("Residue",&residue);
  
  read_str_int("UseSme",&useSme);
  
  if(useSme)
    {
      read_str_int("NApe",&nApe);
      read_str_double("ApeAlpha",&apeAlpha);
      read_str_int("NSme",&nSme);
      read_str_double("KappaSme",&kappaSme);
    }
  
  read_str_int("NGaugeConf",&ngauge_conf);
  
  glb_conf=nissa_malloc("glb_conf",locVol+bord_vol+edge_vol,quad_su3);
  if(useSme)
    ape_conf=nissa_malloc("ape_conf",locVol+bord_vol+edge_vol,quad_su3);
  
  if(cSW!=0.0)
    {
      Cl=nissa_malloc("Cl",locVol,clover_term_t);
      invCl=nissa_malloc("invCl",locVol,inv_clover_term_t);
    }
  
  source_ptr=nissa_malloc("source_ptr",glbSize[0],spincolor*);
  for(int t=0;t<glbSize[0];t++)
    source(t)=nissa_malloc("source",locVol+bord_vol,spincolor);
  
  temp=nissa_malloc("temp",locVol+bord_vol,spincolor);
  
  loc_contr=nissa_malloc("loc_contr",locVol,complex);
  temp_contr=nissa_malloc("temp_contr",glbSize[0],complex);
  
  prop_ptr=nissa_malloc("prop_ptr",2*glbSize[0],spincolor*);
  for(int t=0;t<glbSize[0];t++)
    for(int r=0;r<2;r++)
      prop(t,r)=nissa_malloc("prop",locVol+bord_vol,spincolor);
  
  if(useNewGenerator)
    field_rng_stream.init(seed);
  else
    nissa::start_loc_rnd_gen(seed);
  
  runMode=SEARCHING_CONF;
}

//close the program
void close()
{
  close_input();
  
  master_printf("\n");
  master_printf("Nanalyzed confs: %d\n",nanalyzed_conf);
  master_printf("Total time: %lg s\n",take_time()-init_time);
  if(nanalyzed_conf)
    {
      master_printf("Time per conf: %lg s\n",(take_time()-init_time)/nanalyzed_conf);
      master_printf(" Time to prepare: %lg s, %lg %c\n",prepare_time,prepare_time/(take_time()-init_time)*100,'%');
      master_printf(" Time to prop: %lg s, %lg %c\n",prop_time,prop_time/(take_time()-init_time)*100,'%');
      master_printf(" Time to S0P5: %lg s, %lg %c\n",S0P5_time,S0P5_time/(take_time()-init_time)*100,'%');
      master_printf(" Time to P5P5: %lg s, %lg %c\n",P5P5_time,P5P5_time/(take_time()-init_time)*100,'%');
      master_printf(" Time to disco: %lg s, %lg %c\n",disco_time,disco_time/(take_time()-init_time)*100,'%');
      master_printf(" Time to P5P5_with_ins: %lg s, %lg %c\n",P5P5_with_ins_time,P5P5_with_ins_time/(take_time()-init_time)*100,'%');
    }
  master_printf("\n");
  
  if(cSW!=0.0)
    {
      nissa_free(Cl);
      nissa_free(invCl);
    }
  
  nissa_free(glb_conf);
  if(useSme)
    nissa_free(ape_conf);
  
  for(int t=0;t<glbSize[0];t++)
    for(int r=0;r<2;r++)
      nissa_free(prop(t,r));
  nissa_free(prop_ptr);
  for(int t=0;t<glbSize[0];t++)
    nissa_free(source(t));
  nissa_free(source_ptr);
  nissa_free(temp);
  nissa_free(temp_contr);
  nissa_free(loc_contr);
}

//check if asked to stop or restart
bool asked_to_stop_or_restart()
{
  const bool asked_stop=file_exists("stop");
  master_printf("Asked to stop: %d\n",asked_stop);
  const int asked_restart=file_exists("restart");
  master_printf("Asked to restart: %d\n",asked_restart);
  
  return asked_stop or asked_restart;
}

//check that enough time is available
bool enough_time()
{
  const double passed_time=take_time()-init_time;
  const double remaining_time=wall_time-passed_time;
  const double time_per_conf=passed_time/(nanalyzed_conf+1e-10);
  const double time_per_conf_with_tol=time_per_conf*1.1;
  const bool no_analyzed_conf=(nanalyzed_conf==0);
  if(not no_analyzed_conf)
    {
      master_printf("Time per conf: %lg s (%d confs)\n",time_per_conf,nanalyzed_conf);
      master_printf("Time per conf with tol: %lg s\n",time_per_conf_with_tol);
      master_printf("Remaining time: %lg s\n",remaining_time);
    }
  
  const bool enough_remaining_time=no_analyzed_conf or (remaining_time>time_per_conf_with_tol);
  if(no_analyzed_conf)
    master_printf("No configuration analyzed yet, proceeding in any case\n");
  else
    if(enough_remaining_time)
      master_printf("Time is enough, going on\n");
    else
      master_printf("Not enough time\n");
  
  return enough_remaining_time;
}

//check if all confs has been analyzed
bool finished_confs()
{
  return iconf>=ngauge_conf;
}

bool read_conf_path_and_check_outpath_not_exists()
{
  read_str(conf_path,1024);
  
  read_str(outfolder,1024);
  
  master_printf("Considering configuration \"%s\" with output path \"%s\".\n",conf_path,outfolder);
  
  const bool has_to_be_created=not dir_exists(outfolder);
  if(has_to_be_created)
      master_printf(" Output path \"%s\" not present.\n",outfolder);
  else
    master_printf(" Output path \"%s\" already present.\n",outfolder);
  
  return has_to_be_created;
}

bool create_outpath()
{
  const bool created=not create_dir(outfolder);
  
  if(created)
    master_printf(" Output path \"%s\" for conf \"%s\" created\n",outfolder,conf_path);
  else
    master_printf(" Failed to create the output \"%s\" for conf \"%s\".\n",outfolder,conf_path);
  
  return created;
}

bool create_run_file()
{
  create_outpath();

  if(snprintf(run_file,1024,"%s/running",outfolder)<0) crash("writing %s",run_file);
  
  return lock_file.try_lock(run_file);
}

bool read_conf()
{
  read_ildg_gauge_conf(glb_conf,conf_path);
  master_printf("plaq: %+16.16g\n",global_plaquette_lx_conf(glb_conf));
  
  return true;
}

bool check_lock_file()
{
  const bool lock_file_valid=lock_file.check_lock();
  
  if(not lock_file_valid)
    master_printf("Somebody acquired the lock on %s\n",run_file);
  
  return lock_file_valid;
}

void setup_conf()
{
  prepare_time-=take_time();
  if(useSme)
    ape_smear_conf(ape_conf,glb_conf,apeAlpha,nApe,all_other_spat_dirs[0],1);
  
  momentum_t old_theta={0,0,0,0};
  momentum_t theta={1,0,0,0};
  adapt_theta(glb_conf,old_theta,theta,0,0);
  prepare_time+=take_time();
}

void skipConf()
{
  iconf++;
  
  if(useNewGenerator)
    field_rng_stream.skipDrawers<spincolor>(nhits*glbSize[0]);
  else
    for(int i=0;i<nhits*glbSize[0];i++)
      generate_undiluted_source(source(0),RND_Z4,0);
}

void searchConf()
{
  prepare_time-=take_time();
  if(read_conf_path_and_check_outpath_not_exists() and
     create_run_file() and
     read_conf() and
     check_lock_file())
    {
      setup_conf();
      runMode=ANALYZING_CONF;
    }
  else
    skipConf();
  
  if(finished_confs())
    {
      master_printf("Analyzed all confs, exiting\n\n");
      file_touch(stop_path);
      runMode=STOP;
    }
  prepare_time+=take_time();
}

//mark a conf as finished
void mark_finished()
{
  char fin_file[1024];
  if(snprintf(fin_file,1024,"%s/finished",outfolder)<0) crash("writing %s",fin_file);
  file_touch(fin_file);
  
  iconf++;
  nanalyzed_conf++;
  
  runMode=SEARCHING_CONF;
}

void fill_source(const int glbT)
{
  // double tFrT[1];
  // field_rng_stream.drawScalar(tFrT);
  // const int glbT=tFrT[0]*glb_size[0];
  master_printf("Source position: %d\n",glbT);
  
  auto source_filler=field_rng_stream.getDrawer<spincolor>();
  source_filler.fillField(source(glbT));
  
  NISSA_PARALLEL_LOOP(loclx,0,locVol)
    {
      if(glbT==glbCoordOfLoclx[loclx][0])
	for(int id=0;id<NDIRAC;id++)
	  for(int ic=0;ic<NCOL;ic++)
	    z4Transform(source(glbT)[loclx][id][ic]);//BoxMullerTransform(source[loclx][id][ic]);
      else
	spincolor_put_to_zero(source(glbT)[loclx]);
    }
  NISSA_PARALLEL_LOOP_END;
  
  //master_printf("Box-Muller transformation performed\n");
  
  set_borders_invalid(source(glbT));
}

void get_prop(const int& t,const int& r)
{
  spincolor* p=prop(t,r);
  
  safe_dirac_prod_spincolor(temp,(tau3[r]==-1)?Pminus:Pplus,source(t));
  if(cSW==0.0)
    inv_tmD_cg_eoprec(p,NULL,glb_conf,kappa,mass*tau3[r],1000000,residue,temp);
  else
    {
      invert_twisted_clover_term(invCl,mass*tau3[r],kappa,Cl);
      inv_tmclovD_cg_eoprec(p,NULL,glb_conf,kappa,Cl,invCl,cSW,mass*tau3[r],1000000,residue,temp);
    }
  safe_dirac_prod_spincolor(p,(tau3[r]==-1)?Pminus:Pplus,prop(t,r));
}

void compute_conn_contr(complex* conn_contr,int r1,int r2,int glbT,const dirac_matr& gamma)
{
  vector_reset(loc_contr);
  
  NISSA_PARALLEL_LOOP(ivol,0,locVol)
    {
      spincolor gs;
      unsafe_dirac_prod_spincolor(gs,gamma,prop(glbT,r2)[ivol]);
      spincolor_scalar_prod(loc_contr[ivol],prop(glbT,r1)[ivol],gs);
    }
  NISSA_PARALLEL_LOOP_END;
  
  complex glb_contr[glbSize[0]];
  glb_reduce(glb_contr,loc_contr,locVol,glbSize[0],locSize[0],glbCoordOfLoclx[0][0]);
  
  for(int t=0;t<glbSize[0];t++)
    {
      const int dt=(t+glbSize[0]-glbT)%glbSize[0];
      complex_copy(conn_contr[dt],glb_contr[t]);
    }
}

void compute_inserted_contr(double* contr,spincolor* propBw,spincolor* propFw,const dirac_matr& gamma)
{
  NISSA_PARALLEL_LOOP(ivol,0,locVol)
    {
      unsafe_dirac_prod_spincolor(temp[ivol],gamma,propFw[ivol]);
    }
  NISSA_PARALLEL_LOOP_END;
  set_borders_invalid(temp);
  
  complex_vector_glb_scalar_prod(contr,(complex*)propBw,(complex*)temp,locVol*sizeof(spincolor)/sizeof(complex));
}

void analyzeConf()
{
  FILE* disco_contr_file=open_file(combine("%s/disco_contr",outfolder),"w");
  FILE* P5P5_contr_file=open_file(combine("%s/P5P5_contr",outfolder),"w");
  FILE* P5P5_contr_sme_file=nullptr;
  if(useSme)
    P5P5_contr_sme_file=open_file(combine("%s/P5P5_contr_sme",outfolder),"w");
  FILE* S0P5_contr_file=open_file(combine("%s/S0P5_contr",outfolder),"w");
  FILE* P5P5_SS_contr_file=open_file(combine("%s/P5P5_SS_contr",outfolder),"w");
  FILE* P5P5_0S_contr_file=open_file(combine("%s/P5P5_0S_contr",outfolder),"w");
  
  if(cSW!=0.0) clover_term(Cl,cSW,glb_conf);
  
  for(int ihit=0;ihit<nhits;ihit++)
    {
      prepare_time-=take_time();
      
      //Fill sources
      for(int glbT=0;glbT<glbSize[0];glbT++)
	{
	  // Fill the source
	  if(useNewGenerator)
	    fill_source(glbT);
	  else
	    generate_undiluted_source(source(glbT),RND_Z4,glbT);
	}
      prepare_time+=take_time();
      
      for(int isme=0;isme<(useSme?2:1);isme++)
	{
	  //Compute props
	  prop_time-=take_time();
	  for(int glbT=0;glbT<glbSize[0];glbT++)
	    for(int r=0;r<2;r++)
	      get_prop(glbT,r);
	  prop_time+=take_time();
	  
	  //Compute disco
	  disco_time-=take_time();
	  complex disco_contr[2][glbSize[0]];
	  for(int r=0;r<2;r++)
	    for(int glbT=0;glbT<glbSize[0];glbT++)
	      {
		prop_multiply_with_gamma(temp,5,prop(glbT,r),-1);
		complex_vector_glb_scalar_prod(disco_contr[r][glbT],(complex*)source(glbT),(complex*)temp,locVol*sizeof(spincolor)/sizeof(complex));
	      }
	  disco_time+=take_time();
	  
	  //Print disco
	  for(int r=0;r<2;r++)
	    {
	      master_fprintf(disco_contr_file,"\n# sme %d , hit %d , r %d\n\n",isme,ihit,r);
	      for(int t=0;t<glbSize[0];t++)
		master_fprintf(disco_contr_file,"%.16lg %.16lg\n",disco_contr[r][t][RE],disco_contr[r][t][IM]);
	    }
	  
	  //Compute S0P5
	  S0P5_time-=take_time();
	  complex S0P5_contr[2][2][glbSize[0]];
	  for(int r1=0;r1<2;r1++)
	    for(int r2=0;r2<2;r2++)
	      for(int glbT=0;glbT<glbSize[0];glbT++)
		{
		  compute_conn_contr(temp_contr,r1,r2,glbT,base_gamma[5]);
		  
		  complex_put_to_zero(S0P5_contr[r1][r2][glbT]);
		  for(int t=0;t<glbSize[0];t++)
		    complex_summassign(S0P5_contr[r1][r2][glbT],temp_contr[t]);
		}
	  S0P5_time+=take_time();
	  
	  //Print S0P5
	  for(int r1=0;r1<2;r1++)
	    for(int r2=0;r2<2;r2++)
	      {
		master_fprintf(S0P5_contr_file,"\n# sme %d , hit %d , r1 %d , r2 %d\n\n",isme,ihit,r1,r2);
		for(int t=0;t<glbSize[0];t++)
		  master_fprintf(S0P5_contr_file,"%.16lg %.16lg\n",S0P5_contr[r1][r2][t][RE],S0P5_contr[r1][r2][t][IM]);
	      }
	  
	  //Compute P5P5_SS and 0S
	  P5P5_with_ins_time-=take_time();
	  complex P5P5_SS_contr[2][2][2][2][glbSize[0]];
	  memset(P5P5_SS_contr,0,sizeof(complex)*2*2*2*2*glbSize[0]);
	  complex P5P5_0S_contr[2][2][2][glbSize[0]];
	  memset(P5P5_0S_contr,0,sizeof(complex)*2*2*2*glbSize[0]);
	  for(int dT=0;dT<glbSize[0];dT++)
	    for(int glbT1=0;glbT1<glbSize[0];glbT1++)
	      {
		int glbT2=(glbT1+dT)%glbSize[0];
		complex c[2][2];
		for(int r1=0;r1<2;r1++)
		  for(int r2=0;r2<2;r2++)
		    compute_inserted_contr(c[r1][r2],prop(glbT2,!r2),prop(glbT1,r1),base_gamma[5]);
		
		//prop with other source, to compute 0S
		complex spect[2];
		for(int r=0;r<2;r++)
		  compute_inserted_contr(spect[r],source(glbT2),prop(glbT1,r),base_gamma[5]);
		
		//SS
		for(int r1f=0;r1f<2;r1f++)
		  for(int r2f=0;r2f<2;r2f++)
		    for(int r1b=0;r1b<2;r1b++)
		      for(int r2b=0;r2b<2;r2b++)
			complex_summ_the_conj1_prod(P5P5_SS_contr[r1b][r2b][r1f][r2f][dT],c[r1b][r2b],c[r1f][r2f]);
		
		//0S
		for(int r1f=0;r1f<2;r1f++)
		  for(int r2f=0;r2f<2;r2f++)
		    for(int rb=0;rb<2;rb++)
		      complex_summ_the_conj1_prod(P5P5_0S_contr[rb][r1f][r2f][dT],spect[rb],c[r1f][r2f]);
	      }
	  
	  //Print P5P5_SS
	  for(int r1b=0;r1b<2;r1b++)
	    for(int r2b=0;r2b<2;r2b++)
	      for(int r1f=0;r1f<2;r1f++)
		for(int r2f=0;r2f<2;r2f++)
		  {
		    master_fprintf(P5P5_SS_contr_file,"\n# sme %d , hit %d , r1b %d , r2b %d , r1f %d , r2f %d\n\n",isme,ihit,r1b,r2b,r1f,r2f);
		    for(int t=0;t<glbSize[0];t++)
		      master_fprintf(P5P5_SS_contr_file,"%.16lg %.16lg\n",
				     P5P5_SS_contr[r1b][r2b][r1f][r2f][t][RE]/glbSize[0],
				     P5P5_SS_contr[r1b][r2b][r1f][r2f][t][IM]/glbSize[0]);
		  }
	  
	  //Print P5P5_0S
	  for(int rb=0;rb<2;rb++)
	    for(int r1f=0;r1f<2;r1f++)
	      for(int r2f=0;r2f<2;r2f++)
		{
		  master_fprintf(P5P5_0S_contr_file,"\n# sme %d , hit %d , rb %d , r1f %d , r2f %d\n\n",isme,ihit,rb,r1f,r2f);
		  for(int t=0;t<glbSize[0];t++)
		    master_fprintf(P5P5_0S_contr_file,"%.16lg %.16lg\n",
				   P5P5_0S_contr[rb][r1f][r2f][t][RE]/glbSize[0],
				   P5P5_0S_contr[rb][r1f][r2f][t][IM]/glbSize[0]);
		}
	  
	  P5P5_with_ins_time+=take_time();
	  
	  //Compute and print P5P5
	  P5P5_time-=take_time();
	  for(int withoutWithSme=0;withoutWithSme<(useSme?2:1);withoutWithSme++)
	    {
	      if(withoutWithSme)
		for(int r1=0;r1<2;r1++)
		  for(int glbT=0;glbT<glbSize[0];glbT++)
		    gaussian_smearing(prop(glbT,r1),prop(glbT,r1),ape_conf,kappaSme,nSme);
	      
	      for(int r1=0;r1<2;r1++)
		for(int r2=0;r2<2;r2++)
		  for(int glbT=0;glbT<glbSize[0];glbT++)
		    {
		      compute_conn_contr(temp_contr,r1,r2,glbT,base_gamma[0]);
		      
		      FILE* outFile=withoutWithSme?P5P5_contr_sme_file:P5P5_contr_file;
		      
		      master_fprintf(outFile,"\n# sme %d , hit %d , r1 %d , r2 %d\n\n",isme,ihit,r1,r2);
		      for(int t=0;t<glbSize[0];t++)
			master_fprintf(outFile,"%.16lg %.16lg\n",temp_contr[t][RE],temp_contr[t][IM]);
		    }
	    }
	  P5P5_time+=take_time();
	}
      
      // Smear it
      if(useSme)
	for(int glbT=0;glbT<glbSize[0];glbT++)
	  gaussian_smearing(source(glbT),source(glbT),ape_conf,kappaSme,nSme);
    }
  
  close_file(disco_contr_file);
  close_file(P5P5_contr_file);
  if(useSme)
    close_file(P5P5_contr_sme_file);
  close_file(S0P5_contr_file);
  close_file(P5P5_SS_contr_file);
  close_file(P5P5_0S_contr_file);
  
  mark_finished();
  
  if(asked_to_stop_or_restart()
     or not enough_time())
    runMode=STOP;
  else
    runMode=SEARCHING_CONF;
}

void in_main(int narg,char **arg)
{
  do
    switch(runMode)
      {
      case INITIALIZING:
	init_simulation(narg,arg);
	break;
      case SEARCHING_CONF:
	searchConf();
	break;
      case ANALYZING_CONF:
	analyzeConf();
	break;
      case STOP:
	break;
      }
  while(runMode!=STOP);
  
  //close the simulation
  close();
}

int main(int narg,char **arg)
{
  init_nissa_threaded(narg,arg,in_main);
  close_nissa();
  
  return 0;
}
