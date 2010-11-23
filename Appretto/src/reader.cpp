#pragma once

#include <lemon.h>
#include <lime.h>
#include "endianess.cpp"

using namespace std;

//Write a vector of doubles
void read_double_vector(LemonReader *reader,void *data,int ndoubles_per_site)
{
  //take initial time
  double tic;
  if(debug)
    {
      MPI_Barrier(cart_comm);
      tic=MPI_Wtime();
    }

  int loc_ndoubles_tot=ndoubles_per_site*loc_vol;

  char *header_type=NULL;
  int glb_dims[4]={glb_size[0],glb_size[1],glb_size[2],glb_size[3]};
  int scidac_mapping[4]={0,1,2,3};

  //swap the endianess if needed
  double *swapped_data;
  if(big_endian) swapped_data=new double[loc_ndoubles_tot];
  else swapped_data=(double*) data;

  while(lemonReaderNextRecord(reader)!=LIME_EOF)
    {
      header_type=lemonReaderType(reader);
      if(rank==0) cout<<"found record: "<<header_type<<endl;
      if (strcmp("scidac-binary-data",header_type)==0)
	{
	  int bytes=lemonReaderBytes(reader);
	  int bytes_supposed=ndoubles_per_site*sizeof(double)*glb_vol;
	  if(bytes!=bytes_supposed)
	    {
	      cerr<<"Opsss! The record contain "<<bytes<<" bytes, and it is supposed to contain: "<<bytes_supposed<<endl;
	      MPI_Abort(MPI_COMM_WORLD,1);
	      MPI_Finalize();
	    }
	  lemonReadLatticeParallelMapped(reader,swapped_data,ndoubles_per_site*sizeof(double),glb_dims,scidac_mapping);
	  //swap the endianess
	  if(big_endian) revert_endianess_double_vector((double*)data,swapped_data,loc_ndoubles_tot);
	}
    }

  if(rank==0) cout<<"Data read!"<<endl;

  if(big_endian) delete[] swapped_data;
}

//Read a whole spincolor
void read_spincolor(char *path,spincolor *spinor)
{
  //Open the file
  MPI_File *reader_file=new MPI_File;
  MPI_File_open(cart_comm,path,MPI_MODE_RDONLY,MPI_INFO_NULL,reader_file);
  LemonReader *reader=lemonCreateReader(reader_file,cart_comm);

  read_double_vector(reader,spinor,sizeof(spincolor)/8);

  lemonDestroyReader(reader);
  MPI_File_close(reader_file);
}
