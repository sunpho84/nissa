#pragma once

//create an nissa reader
nissa_reader *nissa_reader_create()
{
  nissa_reader *reader=nissa_malloc("reader",1,nissa_reader);
  
  reader->open=0;
  reader->reading=0;
  reader->buf=NULL;
  reader->lemon_reader=NULL;
  reader->reader_file=nissa_malloc("mpi_file",1,MPI_File);
  return reader;
}

//close an nissa reader
void nissa_reader_close(nissa_reader *reader)
{
  if(reader->open) 
    {
      lemonDestroyReader(reader->lemon_reader);
      MPI_File_close(reader->reader_file);
      reader->open=0;
    }
}

//destroy an nissa reader
void nissa_reader_destroy(nissa_reader *reader)
{
  if(reader->open) nissa_reader_close(reader);
  nissa_free(reader->reader_file);
  if(reader->buf!=NULL) nissa_free(reader->buf);
  if(reader->open) nissa_reader_close(reader);

  nissa_free(reader);
}

//open a file through an nissa reader
void nissa_reader_open(nissa_reader *reader,char *path)
{
  if(reader->open) nissa_reader_close(reader);
  
  if(MPI_File_open(cart_comm,path,MPI_MODE_RDONLY,MPI_INFO_NULL,reader->reader_file)!=MPI_SUCCESS)
    crash("Couldn't open for reading the file: '%s'",path);
  reader->open=1;

  reader->lemon_reader=lemonCreateReader(reader->reader_file,cart_comm);
}

//search a particular record in a file
int nissa_reader_search_record(nissa_reader *reader,const char *expected_record)
{
  int found=0;
  while(found==0 && lemonReaderNextRecord(reader->lemon_reader)!=LEMON_EOF)
    {
      char *header=lemonReaderType(reader->lemon_reader);

      verbosity_lv3_master_printf("found record: %s\n",header);
      if(strcmp(expected_record,header)==0) found=1;
    }
  
  return found;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//start reading current record of a file
//actually it reads and proceeds as it is not clear that non blocking function works really
void nissa_reader_start_reading_current_record(void *out,nissa_reader *reader,int max_nbytes_per_site)
{
  if(!(reader->open)) crash("Ops! Trying to read from a not open reader!");
  
  uint64_t nbytes=lemonReaderBytes(reader->lemon_reader);
  int nbytes_per_site=nbytes/glb_vol;
  if(nbytes_per_site>max_nbytes_per_site)
    crash("Opsss! The file contain %d bytes per site and it is supposed to contain not more than %d!",
	  nbytes_per_site,max_nbytes_per_site);
  
  int glb_dims[4]={glb_size[0],glb_size[3],glb_size[2],glb_size[1]};
  int scidac_mapping[4]={0,3,2,1};
  if(lemonReadLatticeParallelMapped(reader->lemon_reader,out,nbytes_per_site,glb_dims,scidac_mapping)!=
     LEMON_SUCCESS) crash("Error starting to read from a file!");
  reader->reading=1;
  reader->nbytes_per_site=nbytes_per_site;
}

//create the reader, open the file, search the record and start reading
nissa_reader *nissa_reader_start_reading(void *out,char *filename,const char *record_name,int max_bytes_per_site)
{
  nissa_reader *reader=nissa_reader_create();
  nissa_reader_open(reader,filename);
  int found=nissa_reader_search_record(reader,record_name);
  if(!found) crash("Error, record %s not found.\n",record_name);
  nissa_reader_start_reading_current_record(out,reader,max_bytes_per_site);
  
  return reader;
}

//start reading a real vector
nissa_reader *start_reading_real_vector(double *out,char *path,const char *expected_record,int nreals_per_site)
{return nissa_reader_start_reading(out,path,expected_record,nreals_per_site*sizeof(double));}

//start reading a color
nissa_reader *start_reading_color(color *out,char *path)
{return start_reading_real_vector((double*)out,path,"scidac-binary-data",nreals_per_color);}

//start reading a spincolor
nissa_reader *start_reading_spincolor(spincolor *out,char *path)
{return start_reading_real_vector((double*)out,path,"scidac-binary-data",nreals_per_spincolor);}

//start reading a colorspinspin
nissa_reader **start_reading_colorspinspin(colorspinspin *out,char *base_path,char *end_path)
{
  nissa_reader **readers=(nissa_reader**)malloc(4*sizeof(nissa_reader*));
  for(int so=0;so<4;so++)
    {
      char filename[1024];
      if(end_path!=NULL) sprintf(filename,"%s.0%d.%s",base_path,so,end_path);
      else sprintf(filename,"%s.0%d",base_path,so);
      readers[so]=nissa_reader_start_reading(out,filename,"scidac-binary-data",nreals_per_spincolor*sizeof(double));
    }
  
  return readers;
}

//start reading a gauge conf
nissa_reader *start_reading_gauge_conf(quad_su3 *out,char *path)
{return start_reading_real_vector((double*)out,path,"ildg-binary-data",nreals_per_quad_su3);}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

//finalize (wait) to read a certain record
//actually it doesn't wait because reading is always done nonblocky
void nissa_reader_finalize_reading_current_record(nissa_reader *reader)
{
  if(!reader->reading) crash("Error! Waiting to finish reading from a reader that is not reading!");
  
  //if(lemonFinishReading(reader->lemon_reader)!=LEMON_SUCCESS) crash("Error waiting to finish reading!");
  reader->reading=0;
}

//read the checksum
void read_checksum(checksum check_read,nissa_reader *reader)
{
  int found=nissa_reader_search_record(reader,"scidac-checksum");
  if(found)
    {
      uint64_t nbytes=lemonReaderBytes(reader->lemon_reader);
      char *mess=(char*)calloc(nbytes+1,sizeof(char));
      if(lemonReaderReadData(mess,&nbytes,reader->lemon_reader)!=LEMON_SUCCESS) crash("Error while reading checksum");
      sscanf(strstr(mess,"<suma>")+6,"%x",&(check_read[0]));
      sscanf(strstr(mess,"<sumb>")+6,"%x",&(check_read[1]));
      
      free(mess);
    }
  else check_read[0]=check_read[1]=0;
}

//wait to finish reading, read the checksum and destroy the reader
void nissa_reader_finalize_reading(checksum read_check,nissa_reader *reader)
{
  nissa_reader_finalize_reading_current_record(reader);
  read_checksum(read_check,reader);  
  nissa_reader_close(reader);
  nissa_reader_destroy(reader);
}

//finalization of a real vector: change endianess if needed
void finalize_reading_real_vector(double *out,nissa_reader *reader,int nreals_per_site)
{
  int nbytes_per_site_float=nreals_per_site*sizeof(float);
  int nbytes_per_site_double=nreals_per_site*sizeof(double);
  int nbytes_per_site_read=reader->nbytes_per_site;
  
  //read the checksum
  checksum read_check={0,0};
  nissa_reader_finalize_reading(read_check,reader);
  if(read_check[0]!=0||read_check[1]!=0)
    {
      master_printf("Checksums read:      %#010x %#010x\n",read_check[0],read_check[1]);

      //compute checksum
      checksum comp_check;
      checksum_compute_ildg_data(comp_check,out,nbytes_per_site_read);
      
      //print the comparison between checksums
      master_printf("Checksums computed:  %#010x %#010x\n",comp_check[0],comp_check[1]);
      if((read_check[0]!=comp_check[0])||(read_check[1]!=comp_check[1])) master_printf("Warning, checksums do not agree!\n");
    }
  else master_printf("Data checksum not found.\n");
    
  if(nbytes_per_site_read!=nbytes_per_site_float && nbytes_per_site_read!=nbytes_per_site_double)
    crash("Opsss! The file contain %d bytes per site and it is supposed to contain: %d (single) or %d (double)",
	  nbytes_per_site_read,nbytes_per_site_float,nbytes_per_site_double);
  
  int loc_nreals_tot=nreals_per_site*loc_vol;
  
  if(nbytes_per_site_read==nbytes_per_site_float) //cast to double changing endianess if needed
    if(little_endian) floats_to_doubles_changing_endianess((double*)out,(float*)out,loc_nreals_tot);
    else floats_to_doubles_same_endianess((double*)out,(float*)out,loc_nreals_tot);
  else //swap the endianess if needed
    if(little_endian) doubles_to_doubles_changing_endianess((double*)out,(double*)out,loc_nreals_tot);
}  

//reorder a read color
void reorder_read_color(color *c)
{
  int *order=nissa_malloc("order",loc_vol,int);
  
  int x[4];
  for(x[0]=0;x[0]<loc_size[0];x[0]++)
    for(x[1]=0;x[1]<loc_size[1];x[1]++)
      for(x[2]=0;x[2]<loc_size[2];x[2]++)
	for(x[3]=0;x[3]<loc_size[3];x[3]++)
	  {
	    int isour=x[1]+loc_size[1]*(x[2]+loc_size[2]*(x[3]+loc_size[3]*x[0]));
	    int idest=loclx_of_coord(x);
	    order[isour]=idest;
	  }

  reorder_vector((char*)c,order,loc_vol,sizeof(color));
  nissa_free(order);
}

//reorder a read spincolor
void reorder_read_spincolor(spincolor *sc)
{
  int *order=nissa_malloc("order",loc_vol,int);
  
  int x[4];
  for(x[0]=0;x[0]<loc_size[0];x[0]++)
    for(x[1]=0;x[1]<loc_size[1];x[1]++)
      for(x[2]=0;x[2]<loc_size[2];x[2]++)
	for(x[3]=0;x[3]<loc_size[3];x[3]++)
	  {
	    int isour=x[1]+loc_size[1]*(x[2]+loc_size[2]*(x[3]+loc_size[3]*x[0]));
	    int idest=loclx_of_coord(x);
	    order[isour]=idest;
	  }

  reorder_vector((char*)sc,order,loc_vol,sizeof(spincolor));
  nissa_free(order);
}

//reorder a read colorspinspin
void reorder_read_colorspinspin(colorspinspin *css)
{
  int *order=nissa_malloc("order",loc_vol*48,int);
  
  int x[4];
  for(x[0]=0;x[0]<loc_size[0];x[0]++)
    for(x[1]=0;x[1]<loc_size[1];x[1]++)
      for(x[2]=0;x[2]<loc_size[2];x[2]++)
	for(x[3]=0;x[3]<loc_size[3];x[3]++)
	  {
	    int ivsour=x[1]+loc_size[1]*(x[2]+loc_size[2]*(x[3]+loc_size[3]*x[0]));
	    int ivdest=loclx_of_coord(x);
	    for(int so=0;so<4;so++)
	      for(int si=0;si<4;si++)
		for(int c=0;c<3;c++)
		  {
		    int isour=12*(so*loc_vol+ivsour)+3*si+c;
		    int idest=so+4*(si+4*(c+3*ivdest));
		    
		    order[isour]=idest;
		  }
	  }
  
  reorder_vector((char*)css,order,loc_vol*48,sizeof(complex));
  nissa_free(order);
}

//reorder a read gauge conf
void reorder_read_ildg_gauge_conf(quad_su3 *conf)
{
  int *order=nissa_malloc("order",loc_vol*4,int);
  
  int x[4];
  for(x[0]=0;x[0]<loc_size[0];x[0]++)
    for(x[1]=0;x[1]<loc_size[1];x[1]++)
      for(x[2]=0;x[2]<loc_size[2];x[2]++)
	for(x[3]=0;x[3]<loc_size[3];x[3]++)
	  {
	    int isour=x[1]+loc_size[1]*(x[2]+loc_size[2]*(x[3]+loc_size[3]*x[0]));
	    int idest=loclx_of_coord(x);
	  
	    for(int i=0;i<3;i++) order[4*isour+i]=4*idest+i+1;
	    order[4*isour+3]=4*idest;
	  }
  
  reorder_vector((char*)conf,order,4*loc_vol,sizeof(su3));
  nissa_free(order);
}

//finalize reading a color
void finalize_reading_color(color *c,nissa_reader *reader)
{
  finalize_reading_real_vector((double*)c,reader,nreals_per_color);
  reorder_read_color(c);
}

//finalize reading a spincolor
void finalize_reading_spincolor(spincolor *sc,nissa_reader *reader)
{
  finalize_reading_real_vector((double*)sc,reader,nreals_per_spincolor);
  reorder_read_spincolor(sc);
}

//finalize reading a colorspinspin
void finalize_reading_colorspinspin(colorspinspin *css,nissa_reader **reader)
{
  for(int i=0;i<4;i++) finalize_reading_real_vector((double*)css,reader[i],nreals_per_spincolor);
  reorder_read_colorspinspin(css);
}

//finalize reading the reading a gauge conf
void finalize_reading_gauge_conf(quad_su3 *conf,nissa_reader *reader)
{
  finalize_reading_real_vector((double*)conf,reader,nreals_per_quad_su3);
  reorder_read_ildg_gauge_conf(conf);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////

//read a binary blob
int read_binary_blob(void *out,char *path,const char *expected_record,int nmax_bytes_per_site)
{
  //Take inital time
  double start_time=take_time();
  
  nissa_reader *reader=nissa_reader_start_reading(out,path,expected_record,nmax_bytes_per_site);
  checksum check;
  nissa_reader_finalize_reading(check,reader);
  int nbytes_per_site=reader->nbytes_per_site;
  
  verbosity_lv2_master_printf("Total time elapsed in reading %Ld bytes: %f s\n",nbytes_per_site*(long long int)glb_vol,take_time()-start_time);
  
  return nbytes_per_site;
}

//read a real vector
void read_real_vector(double *out,char *path,const char *expected_record,int nreals_per_site)
{
  master_printf("Reading vector: %s\n",path);
  
  //Take inital time
  double start_time=take_time();
  nissa_reader *reader=start_reading_real_vector(out,path,expected_record,nreals_per_site);  
  
  finalize_reading_real_vector(out,reader,nreals_per_site);
  
  verbosity_lv2_master_printf("Total time elapsed including possible conversion: %f s\n",take_time()-start_time);
}

//read a color
void read_color(color *c,char *path)
{
  read_real_vector((double*)c,path,"scidac-binary-data",nreals_per_color);
  reorder_read_color(c);
}

//read a spincolor
void read_spincolor(spincolor *sc,char *path)
{
  read_real_vector((double*)sc,path,"scidac-binary-data",nreals_per_spincolor);
  reorder_read_spincolor(sc);
}

//read a colorspinspin
void read_colorspinspin(colorspinspin *css,char *base_path,char *end_path)
{
  spincolor *temp=nissa_malloc("temp",loc_vol,spincolor);
  for(int so=0;so<4;so++)
    {
      char filename[1024];  
      if(end_path!=NULL) sprintf(filename,"%s.0%d.%s",base_path,so,end_path);
      else sprintf(filename,"%s.0%d",base_path,so);
      read_spincolor(temp,filename);
      nissa_loc_vol_loop(ivol)
	put_spincolor_into_colorspinspin(css[ivol],temp[ivol],so);
    }
  set_borders_invalid(css);
  nissa_free(temp);
}

//read an su3spinspin
void read_su3spinspin(su3spinspin *ccss,char *base_path,char *end_path)
{
  spincolor *temp=nissa_malloc("temp",loc_vol,spincolor);
  for(int id=0;id<4;id++)
    for(int ic=0;ic<3;ic++)
      {
	char filename[1024];  
	if(end_path!=NULL) sprintf(filename,"%s.%02d.%s",base_path,id*3+ic,end_path);
	else sprintf(filename,"%s.%02d",base_path,id*3+ic);
	read_spincolor(temp,filename);
	nissa_loc_vol_loop(ivol)
	  put_spincolor_into_su3spinspin(ccss[ivol],temp[ivol],id,ic);
      }
  set_borders_invalid(ccss);
  nissa_free(temp);
}

//read a gauge conf
void read_ildg_gauge_conf(quad_su3 *conf,char *path)
{
  master_printf("\nReading configuration from file: %s\n",path);
  read_real_vector((double*)conf,path,"ildg-binary-data",nreals_per_quad_su3);
  master_printf("Configuration read!\n\n");
  reorder_read_ildg_gauge_conf(conf);
  set_borders_invalid(conf);
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//Read a spincolor and reconstruct the doublet
void read_tm_spincolor_reconstructing(spincolor **out,spincolor *temp,char *path,quad_su3 *conf,double kappa,double mu)
{
  int all=0;
  if(temp==NULL)
    {
      temp=nissa_malloc("temp",loc_vol+bord_vol,spincolor);
      all=1;
    }
  read_spincolor(temp,path);
  set_borders_invalid(temp);

  reconstruct_tm_doublet(out[0],out[1],conf,kappa,mu,temp);

  if(all) nissa_free(temp);
}  

//Read 4 spincolor and reconstruct them
void read_tm_colorspinspin_reconstructing(colorspinspin **css,char *base_path,char *end_path,quad_su3 *conf,double kappa,double mu)
{
  double start_time=take_time();
  
  char filename[1024];
  spincolor *sc[2]={nissa_malloc("sc1",loc_vol,spincolor),nissa_malloc("sc3",loc_vol,spincolor)};
  spincolor *temp=nissa_malloc("temp",loc_vol+bord_vol,spincolor);

  //Read the four spinor
  for(int id_source=0;id_source<4;id_source++) //dirac index of source
    {
      if(end_path!=NULL) sprintf(filename,"%s.0%d.%s",base_path,id_source,end_path);
      else sprintf(filename,"%s.0%d",base_path,id_source);
      read_tm_spincolor_reconstructing(sc,temp,filename,conf,kappa,mu);
      
      //Switch the spincolor into the colorspin. 
      nissa_loc_vol_loop(ivol)
	{
	  put_spincolor_into_colorspinspin(css[0][ivol],sc[0][ivol],id_source);
	  put_spincolor_into_colorspinspin(css[1][ivol],sc[1][ivol],id_source);
	}
    }

  verbosity_lv1_master_printf("Time elapsed in reading file '%s': %f s\n",base_path,take_time()-start_time);

  set_borders_invalid(css[0]);
  set_borders_invalid(css[1]);
  
  //Destroy the temp
  nissa_free(sc[0]);
  nissa_free(sc[1]);
  nissa_free(temp);
}

//read an ildg conf and split it into e/o parts
void read_ildg_gauge_conf_and_split_into_eo_parts(quad_su3 **eo_conf,char *path)
{
  quad_su3 *lx_conf=nissa_malloc("temp_conf",loc_vol,quad_su3);
  read_ildg_gauge_conf(lx_conf,path);
  split_lx_conf_into_eo_parts(eo_conf,lx_conf);
  nissa_free(lx_conf);

  master_printf("plaq: %.18g\n",global_plaquette_eo_conf(eo_conf));
}
