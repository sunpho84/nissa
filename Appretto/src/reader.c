#pragma once

#include <lemon.h>
#include "endianess.c"
#include "communicate.c"
#include "global.c"
#include "dirac_operator.c"
#include "gaugeconf.c"

//Read from the argument path a maximal amount of data nbyes_per_site
//return the real read amount of bytes
int read_binary_blob(char *data_out,char *path,const char *expected_record,int max_nbytes_per_site)
{
  int nbytes_per_site=0;
  
  //Take inital time
  double tic;
  if(debug>1)
    {
      MPI_Barrier(cart_comm);
      tic=MPI_Wtime();
    }

  //Open the file
  MPI_File *reader_file=(MPI_File*)malloc(sizeof(MPI_File));
  int ok=MPI_File_open(cart_comm,path,MPI_MODE_RDONLY,MPI_INFO_NULL,reader_file);
  if(ok!=MPI_SUCCESS && rank==0)
    {
      fprintf(stderr,"Couldn't open for reading the file: '%s'\n",path);
      fflush(stderr);
      MPI_Abort(cart_comm,1);
    }

  LemonReader *reader=lemonCreateReader(reader_file,cart_comm);
  char *header=NULL;

  int read=0;
  while(lemonReaderNextRecord(reader)!=LEMON_EOF && read==0)
    {
      header=lemonReaderType(reader);
      if(rank==0 && debug>1) printf("found record: %s\n",header);
      if(strcmp(expected_record,header)==0)
	{
	  uint64_t nbytes=lemonReaderBytes(reader);
	  nbytes_per_site=nbytes/glb_vol;
	  if(nbytes_per_site>max_nbytes_per_site && rank==0)
	    {
	      fprintf(stderr,"Opsss! The file contain %Ld bytes per site and it is supposed to contain not more than: %d or %d\n",(long long int)nbytes,nbytes_per_site,max_nbytes_per_site);
	      fflush(stderr);
	      MPI_Abort(MPI_COMM_WORLD,1);
	    }
	  
	  //load with timing
	  double tic1;
	  if(debug>1)
	    {
	      MPI_Barrier(cart_comm);
	      tic1=MPI_Wtime();
	    }
	  int glb_dims[4]={glb_size[0],glb_size[3],glb_size[2],glb_size[1]};
	  int scidac_mapping[4]={0,3,2,1};
	  lemonReadLatticeParallelMapped(reader,data_out,nbytes_per_site,glb_dims,scidac_mapping);
	  if(debug>1)
	    {
	      MPI_Barrier(cart_comm);
	      double tac1=MPI_Wtime();
	      if(rank==0) printf("Time elapsed by lemon to read %Ld bytes: %f s\n",(long long int)nbytes,tac1-tic1);
	    }

	  read=1;
	  if(rank==0 && debug>1) printf("Data read!\n");
	}
    }
  
  //Abort if couldn't find th record
  if(read==0 && rank==0)
    {
      fprintf(stderr,"Error: couldn't find binary\n");
      fflush(stderr);
      MPI_Abort(MPI_COMM_WORLD,1);
    }

  lemonDestroyReader(reader);
  MPI_File_close(reader_file);

 if(debug>1)
    {
      MPI_Barrier(cart_comm);
      double tac=MPI_Wtime();

      if(rank==0) printf("Total time elapsed in reading: %f s\n",tac-tic);
    }

 return nbytes_per_site;
}

//Read a vector of nreal_per_site reals number in float or double precision
//change endianess if needed
void read_real_vector(double *out,char *path,const char *header,int nreals_per_site)
{
  //Take inital time
  double tic;
  if(debug>1)
    {
      MPI_Barrier(cart_comm);
      tic=MPI_Wtime();
    }

  int nbytes_per_site_float=nreals_per_site*sizeof(float);
  int nbytes_per_site_double=nreals_per_site*sizeof(double);
  int nbytes_per_site_read=read_binary_blob((char*)out,path,header,nbytes_per_site_double);
  
  if(nbytes_per_site_read!=nbytes_per_site_float && nbytes_per_site_read!=nbytes_per_site_double && rank==0)
    {
      fprintf(stderr,"Opsss! The file %s contain %d bytes per site and it is supposed to contain: %d (single) or %d (double)\n",
	      path,nbytes_per_site_read,nbytes_per_site_float,nbytes_per_site_double);
      fflush(stderr);
      MPI_Abort(MPI_COMM_WORLD,1);
    }
  
  int loc_nreals_tot=nreals_per_site*loc_vol;
  
  if(nbytes_per_site_read==nbytes_per_site_float) //cast to double changing endianess if needed
    if(big_endian) floats_to_doubles_changing_endianess((double*)out,(float*)out,loc_nreals_tot);
    else floats_to_doubles_same_endianess((double*)out,(float*)out,loc_nreals_tot);
  else //swap the endianess if needed
    if(big_endian) doubles_to_doubles_changing_endianess((double*)out,(double*)out,loc_nreals_tot);

 if(debug>1)
    {
      MPI_Barrier(cart_comm);
      double tac=MPI_Wtime();

      if(rank==0) printf("Total time including possible conversion: %f s\n",tac-tic);
    }
}  

//Read a whole spincolor
void read_spincolor(spincolor *out,char *path)
{
  spincolor *temp=(spincolor*)malloc(sizeof(spincolor)*loc_vol);

  read_real_vector((double*)temp,path,"scidac-binary-data",nreals_per_spincolor);

  int x[4],isour,idest;

  for(x[0]=0;x[0]<loc_size[0];x[0]++)
    for(x[1]=0;x[1]<loc_size[1];x[1]++)
      for(x[2]=0;x[2]<loc_size[2];x[2]++)
	for(x[3]=0;x[3]<loc_size[3];x[3]++)
	  {
	    isour=x[1]+loc_size[1]*(x[2]+loc_size[2]*(x[3]+loc_size[3]*x[0]));
	    idest=loclx_of_coord(x);

	    memcpy(out[idest],temp[isour],sizeof(spincolor));
	  }
  
  free(temp);
}  

//Read a spincolor and reconstruct the doublet
void read_spincolor_reconstructing(spincolor **out,spincolor *temp,char *path,quad_su3 *conf,double kappa,double mu)
{
  int all=0;
  if(temp==NULL)
    {
      temp=(spincolor*)malloc(sizeof(spincolor)*(loc_vol+loc_bord));
      all=1;
    }
  read_spincolor(temp,path);
  communicate_lx_spincolor_borders(temp);

  reconstruct_doublet(out[0],out[1],temp,conf,kappa,mu);  

  if(all) free(temp);
}  

//Read 4 spincolor and revert their indexes
void read_colorspinspin(colorspinspin *css,char *base_path,char *end_path)
{
  double tic;
  if(debug)
    {
      MPI_Barrier(cart_comm);
      tic=MPI_Wtime();
    }

  char filename[1024];
  spincolor *sc=(spincolor*)malloc(sizeof(spincolor)*loc_vol);

  //Read the four spinor
  for(int id_source=0;id_source<4;id_source++) //dirac index of source
    {
      if(end_path!=NULL) sprintf(filename,"%s.0%d.%s",base_path,id_source,end_path);
      else sprintf(filename,"%s.0%d",base_path,id_source);
      read_spincolor(sc,filename);
      
      //Switch the spincolor into the colorspin. 
      for(int loc_site=0;loc_site<loc_vol;loc_site++) put_spincolor_into_colorspinspin(css[loc_site],sc[loc_site],id_source);
    }

  if(debug)
    {
      MPI_Barrier(cart_comm);
      double tac=MPI_Wtime();

      if(rank==0) printf("Time elapsed in reading file '%s': %f s\n",base_path,tac-tic);
    }
  
  //Destroy the temp
  free(sc);
}

//Read 4 spincolor and reconstruct them
void read_colorspinspin_reconstructing(colorspinspin **css,char *base_path,char *end_path,quad_su3 *conf,double kappa,double mu)
{
  double tic;
  if(debug)
    {
      MPI_Barrier(cart_comm);
      tic=MPI_Wtime();
    }
  
  char filename[1024];
  spincolor *sc[2]={(spincolor*)malloc(sizeof(spincolor)*loc_vol),(spincolor*)malloc(sizeof(spincolor)*loc_vol)};
  spincolor *temp=(spincolor*)malloc(sizeof(spincolor)*(loc_vol+loc_bord));

  //Read the four spinor
  for(int id_source=0;id_source<4;id_source++) //dirac index of source
    {
      if(end_path!=NULL) sprintf(filename,"%s.0%d.%s",base_path,id_source,end_path);
      else sprintf(filename,"%s.0%d",base_path,id_source);
      read_spincolor_reconstructing(sc,temp,filename,conf,kappa,mu);
      
      //Switch the spincolor into the colorspin. 
      for(int loc_site=0;loc_site<loc_vol;loc_site++)
	{
	  put_spincolor_into_colorspinspin(css[0][loc_site],sc[0][loc_site],id_source);
	  put_spincolor_into_colorspinspin(css[1][loc_site],sc[1][loc_site],id_source);
	}
    }

  if(debug)
    {
      MPI_Barrier(cart_comm);
      double tac=MPI_Wtime();

      if(rank==0) printf("Time elapsed in reading file '%s': %f s\n",base_path,tac-tic);
    }

  //Destroy the temp
  free(sc[0]);
  free(sc[1]);
  free(temp);
}

////////////////////////// gauge configuration loading /////////////////////////////

//Read only the local part of the gauge configuration
void read_local_gauge_conf(quad_su3 *out,char *path)
{
  double tic;
  if(debug)
    {
      MPI_Barrier(cart_comm);
      tic=MPI_Wtime();
    }

  quad_su3 *temp=(quad_su3*)malloc(sizeof(quad_su3)*loc_vol);

  read_real_vector((double*)temp,path,"ildg-binary-data",nreals_per_quad_su3);
  
  int x[4],isour,idest;
  quad_su3 buff;

  for(x[0]=0;x[0]<loc_size[0];x[0]++)
    for(x[1]=0;x[1]<loc_size[1];x[1]++)
      for(x[2]=0;x[2]<loc_size[2];x[2]++)
	for(x[3]=0;x[3]<loc_size[3];x[3]++)
	  {
	    isour=x[1]+loc_size[1]*(x[2]+loc_size[2]*(x[3]+loc_size[3]*x[0]));
	    idest=loclx_of_coord(x);

	    memcpy(buff[0],temp[isour][3],sizeof(su3));
	    memcpy(buff[1],temp[isour][0],sizeof(su3));
	    memcpy(buff[2],temp[isour][1],sizeof(su3));
	    memcpy(buff[3],temp[isour][2],sizeof(su3));
	    
	    //this is to avoid premature overwrite
	    memcpy(out[idest],buff,sizeof(quad_su3));
	  }

  free(temp);
  
  if(debug)
    {
      MPI_Barrier(cart_comm);
      double tac=MPI_Wtime();
      
      if(rank==0) printf("Time elapsed in reading gauge file '%s': %f s\n",path,tac-tic);
    }
}
