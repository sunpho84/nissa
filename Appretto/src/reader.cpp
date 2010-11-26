#pragma once

#include <lemon.h>
#include <lime.h>
#include "endianess.cpp"

using namespace std;

//Read a whole spincolor
void read_spincolor(char *path,spincolor *spinor)
{
  const int nreals_per_site=24;

  //take initial time
  double tic;
  if(debug>1)
    {
      MPI_Barrier(cart_comm);
      tic=MPI_Wtime();
    }

  //Open the file
  MPI_File *reader_file=new MPI_File;
  int ok=MPI_File_open(cart_comm,path,MPI_MODE_RDONLY,MPI_INFO_NULL,reader_file);
  if(ok!=MPI_SUCCESS)
    {
      if(rank==0) cerr<<"Couldn't open for reading the file: '"<<path<<"'"<<endl;
      MPI_Abort(cart_comm,1);
    }

  LemonReader *reader=lemonCreateReader(reader_file,cart_comm);
  char *header_type=NULL;

  bool read=false;
  while(lemonReaderNextRecord(reader)!=LIME_EOF)
    {
      header_type=lemonReaderType(reader);
      if(rank==0 and debug>1) cout<<"found record: "<<header_type<<endl;
      if(strcmp("scidac-binary-data",header_type)==0)
	{
	  int bytes=lemonReaderBytes(reader);
	  int bytes_float=nreals_per_site*sizeof(float)*glb_vol;
	  int bytes_double=nreals_per_site*sizeof(double)*glb_vol;
	  if(bytes!=bytes_float and bytes!=bytes_double)
	    {
	      cerr<<"Opsss! The record contain "<<bytes<<" bytes, and it is supposed to contain: "
		  <<bytes_float<<" or "<<bytes_double<<endl;
	      MPI_Abort(MPI_COMM_WORLD,1);
	      MPI_Finalize();
	    }
	  
	  int loc_nreals_tot=nreals_per_site*loc_vol;
	  
	  int glb_dims[4]={glb_size[0],glb_size[1],glb_size[2],glb_size[3]};
	  int scidac_mapping[4]={0,1,2,3};
	  
	  lemonReadLatticeParallelMapped(reader,spinor,bytes,glb_dims,scidac_mapping);
	  
	  if(bytes==bytes_float) //cast to double changing endianess if needed
	    if(big_endian) floats_to_doubles_changing_endianess((double*)spinor,(float*)spinor,loc_nreals_tot);
	    else floats_to_doubles_same_endianess((double*)spinor,(float*)spinor,loc_nreals_tot);
	  else //swap the endianess if needed
	    if(big_endian) doubles_to_doubles_changing_endianess((double*)spinor,(double*)spinor,loc_nreals_tot);
	  
	  read=true;
	  if(rank==0 and debug>1) cout<<"Data read!"<<endl;
	}
    }
  
  if(read==false)
    {
      if(rank==0) cerr<<"Error: couldn't find binary"<<endl;
      MPI_Abort(MPI_COMM_WORLD,1);
    }

  lemonDestroyReader(reader);
  MPI_File_close(reader_file);

 if(debug>1)
    {
      MPI_Barrier(cart_comm);
      double tac=MPI_Wtime();

      if(rank==0) cout<<"Time elapsed in reading "<<tac-tic<<" s"<<endl;
    }
}

//Read 4 spincolor and revert their indexes
void read_colorspinspin(char *base_path,colorspinspin *css)
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
      sprintf(filename,"%s.0%d",base_path,id_source);
      read_spincolor(filename,sc);
      
      //Switch the spincolor into the colorspin. 
      //In a spinspin the source index runs slower than the sink
      for(int loc_site=0;loc_site<loc_vol;loc_site++)
	for(int ic_sink=0;ic_sink<3;ic_sink++)
	  for(int id_sink=0;id_sink<4;id_sink++) //dirac index of sink
	    {
	      css[loc_site][ic_sink][id_source][id_sink][0]=sc[loc_site][id_sink][ic_sink][0];
	      css[loc_site][ic_sink][id_source][id_sink][1]=sc[loc_site][id_sink][ic_sink][1];
	    }
    }

  if(debug)
    {
      MPI_Barrier(cart_comm);
      double tac=MPI_Wtime();

      if(rank==0) cout<<"Time elapsed in reading file '"<<base_path<<"': "<<tac-tic<<" s"<<endl;
    }
  
  //Destroy the temp
  free(sc);
}
