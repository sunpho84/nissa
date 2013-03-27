#ifdef HAVE_CONFIG_H
 #include "config.h"
#endif

#include "nissa.h"

#include <malloc.h>

#include <firmware/include/personality.h>
#include <spi/include/mu/Descriptor_inlines.h>
#include <spi/include/mu/Addressing_inlines.h>
#include <spi/include/kernel/MU.h>
#include <spi/include/kernel/location.h>

//global barrier for spi
void spi_global_barrier()
{
  if(MUSPI_GIBarrierEnter(&spi_barrier)) crash("while entering spi barrier");
  if(MUSPI_GIBarrierPollWithTimeout(&spi_barrier,60UL*1600000000)) crash("while waiting for barrier to finish");
}

//get the spi coord and grid size
void get_spi_coord()
{
  //get the personality
  Personality_t pers;
  Kernel_GetPersonality(&pers,sizeof(pers));
      
  //copy the coords
  spi_rank_coord[0]=pers.Network_Config.Acoord;
  spi_rank_coord[1]=pers.Network_Config.Bcoord;
  spi_rank_coord[2]=pers.Network_Config.Ccoord;
  spi_rank_coord[3]=pers.Network_Config.Dcoord;
  spi_rank_coord[4]=pers.Network_Config.Ecoord;
}

//set the spi neighboirs using MPI communicators
void set_spi_neighbours()
{
  //loop on the 8 dirs
  for(int mu=0;mu<4;mu++)
    for(int bf=0;bf<2;bf++) //backward(0) or forward(1)
      {
	//send to one dir, receive from the opposite
	coords_5D c;
	MPI_Sendrecv(spi_rank_coord,sizeof(coords_5D),MPI_BYTE,rank_neigh[bf][mu],0,
		     c,sizeof(coords_5D),MPI_BYTE,rank_neigh[!bf][mu],0,
		     cart_comm,MPI_STATUS_IGNORE);
	
	//setup the neighbours
	MUSPI_SetUpDestination(&spi_neigh[!bf][mu],c[0],c[1],c[2],c[3],c[4]);
      }
}

//find coordinates and set neighbours map
void set_spi_geometry()
{
  get_spi_coord();
  set_spi_neighbours();      
}

//initialize the spi communications
void init_spi()
{
  //check not to have initialized
  if(!nissa_spi_inited)
    {
      verbosity_lv2_master_printf("Starting spi\n");
      
      //check that we have parallelized all the dirs
      if(nparal_dir!=4) crash("all directions must be parallelized (only %d resulting parallelized)",nparal_dir);
      
      //check that we do not have more than one process per node
      if(Kernel_ProcessCount()!=1) crash("only one process per node implemented");
      
      //mark as initialized
      nissa_spi_inited=true;
      
      //get coordinates, size and rank in the 5D grid
      set_spi_geometry();
      
      //reset the number of allocated bat
      spi_nallocated_bat=0;
      
      ////////////////////////////////// init the fifos ///////////////////////////////////
      
      //alloc space for the 8 injection fifos
      uint32_t fifo_size=64*8;
      for(int idir=0;idir<8;idir++) spi_fifo[idir]=(uint64_t*)memalign(64,fifo_size);

      //set default attributes for inj fifo
      Kernel_InjFifoAttributes_t fifo_attrs[8];
      for(int idir=0;idir<8;idir++) memset(&fifo_attrs[idir],0,sizeof(Kernel_InjFifoAttributes_t));
      
      //initialize them with default attributes
      uint32_t fifo_id[8]={0,1,2,3,4,5,6,7};
      if(Kernel_AllocateInjFifos(0,&spi_fifo_sg_ptr,8,fifo_id,fifo_attrs)) crash("allocating inj fifos");
      
      //init the MU MMIO for the fifos
      for(int idir=0;idir<8;idir++)
	{
	  //create the memory region
	  Kernel_MemoryRegion_t mem_region;
	  if(Kernel_CreateMemoryRegion(&mem_region,spi_fifo[7-idir],fifo_size)) crash("creating memory region");

	  //initialise the fifos
	  if(Kernel_InjFifoInit(&spi_fifo_sg_ptr,fifo_id[idir],&mem_region,
				(uint64_t)spi_fifo[7-idir]-(uint64_t)mem_region.BaseVa,fifo_size-1)) crash("initializing fifo");
	}
      
      //activate the fifos
      if(Kernel_InjFifoActivate(&spi_fifo_sg_ptr,8,fifo_id,KERNEL_INJ_FIFO_ACTIVATE)) crash("activating fifo");
  
      //init the barrier
      if(MUSPI_GIBarrierInit(&spi_barrier,0)) crash("initializing the barrier");
    }
}

//general set of spi comm
void set_spi_comm(spi_comm_t &in,int buf_size,int *payload_address_offset,int *message_length,MUHWI_Destination *dest,int *rec_payload_offset)
{
  /////////////////////////////////////// allocate buffers ///////////////////////////////////////
  
  in.comm_in_prog=0;
  in.buf_size=buf_size;
  
  in.recv_buf=(char*)memalign(64,buf_size);
  in.send_buf=(char*)memalign(64,buf_size);
  in.descriptors=(MUHWI_Descriptor_t*)memalign(64,8*sizeof(MUHWI_Descriptor_t));
  
  //////////////////////////////// allocate base address table (bat) /////////////////////////////
  
  //allocate the bat entries
  in.bat_id[0]=spi_nallocated_bat+0;
  in.bat_id[1]=spi_nallocated_bat+1;
  spi_nallocated_bat+=2;
  if(Kernel_AllocateBaseAddressTable(0,&in.spi_bat_gr,2,in.bat_id,0)) crash("allocating bat");
  verbosity_lv3_master_printf("number of spi allocated bat: %d\n",spi_nallocated_bat);
  
  //get physical address of receiving buffer
  Kernel_MemoryRegion_t mem_region;
  if(Kernel_CreateMemoryRegion(&mem_region,in.recv_buf,in.buf_size)) crash("creating memory region");
  
  //set the physical address
  if(MUSPI_SetBaseAddress(&in.spi_bat_gr,in.bat_id[0],(uint64_t)in.recv_buf-(uint64_t)mem_region.BaseVa+(uint64_t)mem_region.BasePa))
    crash("setting base address");
  
  //set receive counter bat to MU style atomic PA addr of the receive counter
  if((uint64_t)(&in.recv_counter)&0x7) crash("recv counter not 8 byte aligned");
  if(Kernel_CreateMemoryRegion(&mem_region,(void*)&in.recv_counter,sizeof(uint64_t))) crash("creating memory region");  
  if(MUSPI_SetBaseAddress(&in.spi_bat_gr,in.bat_id[1],
			  MUSPI_GetAtomicAddress((uint64_t)&in.recv_counter-(uint64_t)mem_region.BaseVa+(uint64_t)mem_region.BasePa,
						 MUHWI_ATOMIC_OPCODE_STORE_ADD))) crash("setting base addr");
  
  //reset number of byte to be received
  in.recv_counter=0;

  //////////////////////////////// prepares the descriptors /////////////////////////////
  
  //get the send buffer physical address
  if(Kernel_CreateMemoryRegion(&mem_region,in.send_buf,buf_size)) crash("creating memory region");
  uint64_t send_buf_phys_addr=(uint64_t)in.send_buf-(uint64_t)mem_region.BaseVa+(uint64_t)mem_region.BasePa;
  
  //create one descriptor per direction
  for(int idir=0;idir<8;idir++)
    {
      //reset the info
      MUSPI_Pt2PtDirectPutDescriptorInfo_t dinfo;
      memset(&dinfo,0,sizeof(MUSPI_Pt2PtDirectPutDescriptorInfo_t));
      
      //set the parameters
      dinfo.Base.Payload_Address=send_buf_phys_addr+payload_address_offset[idir];
      dinfo.Base.Message_Length=message_length[idir];
	dinfo.Base.Torus_FIFO_Map=MUHWI_DESCRIPTOR_TORUS_FIFO_MAP_AM|MUHWI_DESCRIPTOR_TORUS_FIFO_MAP_AP|
	  MUHWI_DESCRIPTOR_TORUS_FIFO_MAP_BM|MUHWI_DESCRIPTOR_TORUS_FIFO_MAP_BP|MUHWI_DESCRIPTOR_TORUS_FIFO_MAP_CM|
	  MUHWI_DESCRIPTOR_TORUS_FIFO_MAP_CP|MUHWI_DESCRIPTOR_TORUS_FIFO_MAP_DM|MUHWI_DESCRIPTOR_TORUS_FIFO_MAP_DP|
	  MUHWI_DESCRIPTOR_TORUS_FIFO_MAP_EM|MUHWI_DESCRIPTOR_TORUS_FIFO_MAP_EP;
      dinfo.Base.Dest=dest[idir];
      dinfo.Pt2Pt.Hints_ABCD=0;
      dinfo.Pt2Pt.Misc1=MUHWI_PACKET_USE_DETERMINISTIC_ROUTING|MUHWI_PACKET_DO_NOT_ROUTE_TO_IO_NODE;	
      dinfo.Pt2Pt.Misc2=MUHWI_PACKET_VIRTUAL_CHANNEL_DETERMINISTIC;
      dinfo.Pt2Pt.Skip=8; //for checksumming, skip the header
      dinfo.DirectPut.Rec_Payload_Base_Address_Id=in.bat_id[0];
      dinfo.DirectPut.Rec_Payload_Offset=rec_payload_offset[idir];
      dinfo.DirectPut.Rec_Counter_Base_Address_Id=in.bat_id[1];
      dinfo.DirectPut.Rec_Counter_Offset=0;
      dinfo.DirectPut.Pacing=MUHWI_PACKET_DIRECT_PUT_IS_NOT_PACED;
      
      //create the descriptor
      if(MUSPI_CreatePt2PtDirectPutDescriptor(&in.descriptors[idir],&dinfo)) crash("creating the descriptor");
    }
}

//set up a communicator for lx or eo borders
//first 4 communicate to forward nodes, last four to backward nodes
void set_lx_or_eo_spi_comm(spi_comm_t &in,int nbytes_per_site,int lx_eo)
{
  int div_coeff=(lx_eo==0)?1:2; //dividing coeff
  int tot_buf_size=nbytes_per_site*bord_vol/div_coeff;
  int payload_address_offset[8],message_length[8],rec_payload_offset[8];
  MUHWI_Destination dest[8];
  
  //direction of the halo in receiving node: surface is ordered opposite of halo
  for(int bf=0;bf<2;bf++) 
    for(int mu=0;mu<4;mu++)
      {
	int idir=bf*4+mu;
	
	//set the parameters
	payload_address_offset[idir]=(bord_offset[mu]+bord_vol/2*(!bf))*nbytes_per_site/div_coeff;
	message_length[idir]=bord_dir_vol[mu]*nbytes_per_site/div_coeff;
	dest[idir]=spi_neigh[!bf][mu];
	rec_payload_offset[idir]=(bord_offset[mu]+bord_vol/2*bf)*nbytes_per_site/div_coeff;
      }
  
  set_spi_comm(in,tot_buf_size,payload_address_offset,message_length,dest,rec_payload_offset);
}

void set_lx_spi_comm(spi_comm_t &in,int nbytes_per_site) {set_lx_or_eo_spi_comm(in,nbytes_per_site,0);}
void set_eo_spi_comm(spi_comm_t &in,int nbytes_per_site) {set_lx_or_eo_spi_comm(in,nbytes_per_site,1);}

//wait a communication to finish
void spi_comm_wait(spi_comm_t *in)
{
  GET_THREAD_ID();
  
  if(IS_MASTER_THREAD)
    {
      verbosity_lv3_master_printf("Entering spi comm wait\n");
      
      if(in->comm_in_prog)
	{
	  //wait to send everything
	  while(in->comm_in_prog)
	    {
	      verbosity_lv3_master_printf("Waiting to finish sending data with spi\n");
	      in->comm_in_prog=0;
	      for(int idir=0;idir<8;idir++)
		in->comm_in_prog|=!MUSPI_CheckDescComplete(MUSPI_IdToInjFifo(idir,&spi_fifo_sg_ptr),spi_desc_count[idir]);
	    }
	  
	  spi_global_barrier();
	  
	  //wait to receive everything
	  verbosity_lv3_master_printf("Waiting to finish receiving data with spi\n");
	  while(in->recv_counter>0) verbosity_lv3_master_printf("%d/%d bytes remaining to be received\n",in->recv_counter,in->buf_size);  
	}
      else verbosity_lv3_master_printf("Did not have to wait for any spi comm\n");
    }
  
  //all threads must wait
  thread_barrier(SPI_COMM_WAIT_BARRIER);
}

//start the communications
void spi_start_comm(spi_comm_t *in,int *dir_comm,int tot_size)
{
  GET_THREAD_ID();
  
  //wait for any previous communication to finish and mark as new started 
  spi_comm_wait(in);
  
  if(IS_MASTER_THREAD)
    {
      //reset the counter and wait that all have reset
      spi_global_barrier();
      if(tot_size==-1) in->recv_counter=in->buf_size;
      else in->recv_counter=tot_size;
      spi_global_barrier();
      
      //start the injection
      for(int idir=0;idir<8;idir++)
	if(dir_comm==NULL||dir_comm[idir])
	  {
	    verbosity_lv3_master_printf("Injecting %d\n",idir);
	    spi_desc_count[idir]=MUSPI_InjFifoInject(MUSPI_IdToInjFifo(idir,&spi_fifo_sg_ptr),&in->descriptors[idir]);
	    if(spi_desc_count[idir]>(1ll<<57)) crash("msg_InjFifoInject returned %llu when expecting 1, most likely because there is no room in the fifo",spi_desc_count[idir]);
	  }
      
      in->comm_in_prog=1;
    }
}

//unset everything
void unset_spi_comm(spi_comm_t *in)
{
  GET_THREAD_ID();
  
  //wait for any communication to finish
  spi_comm_wait(in);
  
  //free buffers
  if(IS_MASTER_THREAD)
    {
      free(in->recv_buf);
      free(in->send_buf);
    }
}

/////////////////////////////////////// communicating lx vec ///////////////////////////////////

//fill the sending buf using the data inside an lx vec
void fill_spi_sending_buf_with_lx_vec(spi_comm_t *a,void *vec,int nbytes_per_site)
{
  GET_THREAD_ID();
  
  //check buffer size matching
  if(a->buf_size!=nbytes_per_site*bord_vol) crash("wrong buffer size (%d) for %d large border)",a->buf_size,nbytes_per_site*bord_vol);
  
  //copy one by one the surface of vec inside the sending buffer
  NISSA_PARALLEL_LOOP(ibord,0,bord_vol)
    memcpy(a->send_buf+nbytes_per_site*ibord,(char*)vec+surflx_of_bordlx[ibord]*nbytes_per_site,nbytes_per_site);

  //wait that all threads filled their portion
  thread_barrier(SPI_LX_SENDING_BUF_FILL_BARRIER);
}

//extract the information from receiving buffer and put them inside an lx vec
void fill_lx_bord_with_spi_receiving_buf(void *vec,spi_comm_t *a,int nbytes_per_site)
{
  GET_THREAD_ID();
  
  if(IS_MASTER_THREAD)
    {
      crash_if_borders_not_allocated(vec);
      
      //check buffer size matching
      if(a->buf_size!=nbytes_per_site*bord_vol) crash("wrong buffer size (%d) for %d large border)",a->buf_size,nbytes_per_site*bord_vol);
      
      //the buffer is already ordered as the vec border
      memcpy((char*)vec+loc_vol*nbytes_per_site,a->recv_buf,a->buf_size);
    }
  
  //we do not sync, because typically we will set borders as valid
}

//start communication using an lx border
void spi_start_communicating_lx_borders(int *nrequest,spi_comm_t *a,void *vec,int nbytes_per_site)
{
  if(!check_borders_valid(vec))
    {
      GET_THREAD_ID();
      
      //take time and write some debug output
      if(IS_MASTER_THREAD) tot_nissa_comm_time-=take_time();
      verbosity_lv3_master_printf("Start spi communication of lx borders of %s\n",get_vec_name((void*)vec));
      
      //fill the communicator buffer, start the communication and take time
      fill_spi_sending_buf_with_lx_vec(a,vec,nbytes_per_site);
      spi_start_comm(a);
      if(IS_MASTER_THREAD) tot_nissa_comm_time+=take_time();
      
      (*nrequest)=8;
    }
  else (*nrequest)=0;
}

//finish communicating
void spi_finish_communicating_lx_borders(int *nrequest,void *vec,spi_comm_t *a,int nbytes_per_site)
{
  if(!check_borders_valid(vec) && (*nrequest==8))
    {
      GET_THREAD_ID();
      
      //take time and write some output
      if(IS_MASTER_THREAD) tot_nissa_comm_time-=take_time();
      verbosity_lv3_master_printf("Finish spi communication of lx borders of %s\n",get_vec_name((void*)vec));

      //wait communication to finish, fill back the vector and take time
      spi_comm_wait(a);
      fill_lx_bord_with_spi_receiving_buf(vec,a,nbytes_per_site);
      if(IS_MASTER_THREAD) tot_nissa_comm_time+=take_time();
      
      //set border not valid: this auto sync
      set_borders_valid(vec);
      
      (*nrequest)=0;
    }
}

//merge the two
void spi_communicate_lx_borders(void *vec,spi_comm_t *a,int nbytes_per_site)
{
  if(!check_borders_valid(vec))
    {
      verbosity_lv3_master_printf("Sync spi communication of lx borders of %s\n",get_vec_name((void*)vec));
      
      int nrequest;
      spi_start_communicating_lx_borders(&nrequest,a,vec,nbytes_per_site);
      spi_finish_communicating_lx_borders(&nrequest,vec,a,nbytes_per_site);
    }
}

/////////////////////////////////////// communicating e/o vec //////////////////////////////////////

//fill the sending buf using the data inside an ev or odd vec
void fill_spi_sending_buf_with_ev_or_od_vec(spi_comm_t *a,void *vec,int nbytes_per_site,int eo)
{
  GET_THREAD_ID();
  
  //check buffer size matching
  if(a->buf_size!=nbytes_per_site*bord_volh) crash("wrong buffer size (%d) for %d border)",a->buf_size,nbytes_per_site*bord_volh);
  
  //copy one by one the surface of vec inside the sending buffer
  NISSA_PARALLEL_LOOP(ibord,0,bord_volh)
    memcpy(a->send_buf+ibord*nbytes_per_site,(char*)vec+surfeo_of_bordeo[eo][ibord]*nbytes_per_site,nbytes_per_site);

  //wait that all threads filled their portion
  thread_barrier(SPI_EV_OR_OD_SENDING_BUF_FILL_BARRIER);
}

//extract the information from receiving buffer and put them inside an even or odd vec
void fill_ev_or_od_bord_with_spi_receiving_buf(void *vec,spi_comm_t *a,int nbytes_per_site)
{
  GET_THREAD_ID();
  
  if(IS_MASTER_THREAD)
    {
      crash_if_borders_not_allocated(vec);
      
      //check buffer size matching
      if(a->buf_size!=nbytes_per_site*bord_volh)
	crash("wrong buffer size (%d) for %d border)",a->buf_size,nbytes_per_site*bord_volh);
      
      //the buffer is already ordered as the vec border
      memcpy((char*)vec+loc_volh*nbytes_per_site,a->recv_buf,a->buf_size);
    }
  
  //we do not sync, because typically we will set borders as valid
}

//start communication using an ev or od border
void spi_start_communicating_ev_or_od_borders(int *nrequest,spi_comm_t *a,void *vec,int nbytes_per_site,int eo)
{
  if(!check_borders_valid(vec))
    {
      GET_THREAD_ID();
      
      //take time and output debugging info
      if(IS_MASTER_THREAD) tot_nissa_comm_time-=take_time();
      verbosity_lv3_master_printf("Starting spi communication of ev or od borders of %s\n",get_vec_name((void*)vec));

      //fill the communicator buffer, start the communication and take time
      fill_spi_sending_buf_with_ev_or_od_vec(a,vec,nbytes_per_site,eo);
      spi_start_comm(a);
      if(IS_MASTER_THREAD) tot_nissa_comm_time+=take_time();
      
      (*nrequest)=8;
    }
  else (*nrequest)=0;
}

//finish communicating
void spi_finish_communicating_ev_or_od_borders(int *nrequest,void *vec,spi_comm_t *a,int nbytes_per_site)
{
  if(!check_borders_valid(vec) && (*nrequest==8))
    {
      GET_THREAD_ID();
      
      //take time and make some output
      if(IS_MASTER_THREAD) tot_nissa_comm_time-=take_time();
      verbosity_lv3_master_printf("Finish spi communication of ev or od borders of %s\n",get_vec_name((void*)vec));
	  
      //wait communication to finish, fill back the vector and take time
      spi_comm_wait(a);
      fill_ev_or_od_bord_with_spi_receiving_buf(vec,a,nbytes_per_site);
      if(IS_MASTER_THREAD) tot_nissa_comm_time+=take_time();
      
      set_borders_valid(vec);
      
      (*nrequest)=0;
    }
}

//merge the two
void spi_communicate_ev_or_od_borders(void *vec,spi_comm_t *a,int nbytes_per_site,int eo)
{
  if(!check_borders_valid(vec))
    {
      verbosity_lv3_master_printf("Sync spi communication of ev or od borders of %s\n",get_vec_name((void*)vec));
      
      int nrequest;
      spi_start_communicating_ev_or_od_borders(&nrequest,a,vec,nbytes_per_site,eo);
      spi_finish_communicating_ev_or_od_borders(&nrequest,vec,a,nbytes_per_site);
    }  
}

/////////////////////////////////////// communicating e&o vec //////////////////////////////////////

//fill the sending buf using the data inside an ev and odd vec, using lx style inside buf
void fill_spi_sending_buf_with_ev_and_od_vec(spi_comm_t *a,void **vec,int nbytes_per_site)
{
  GET_THREAD_ID();
  
  //check buffer size matching
  if(a->buf_size!=nbytes_per_site*bord_vol) crash("wrong buffer size (%d) for %d border)",a->buf_size,nbytes_per_site*bord_vol);
  
  //copy one by one the surface of vec inside the sending buffer
  NISSA_PARALLEL_LOOP(ibord_lx,0,bord_vol)
    {
      //convert lx indexing to eo
      int source_lx=surflx_of_bordlx[ibord_lx];
      int par=loclx_parity[source_lx];
      int source_eo=loceo_of_loclx[source_lx];
      memcpy(a->send_buf+ibord_lx*nbytes_per_site,(char*)(vec[par])+source_eo*nbytes_per_site,nbytes_per_site);
    }
  
  //wait that all threads filled their portion
  thread_barrier(SPI_EV_AND_OD_SENDING_BUF_FILL_BARRIER);
}

//extract the information from receiving buffer and put them inside an even or odd vec
void fill_ev_and_od_bord_with_spi_receiving_buf(void **vec,spi_comm_t *a,int nbytes_per_site)
{
  GET_THREAD_ID();
  
  crash_if_borders_not_allocated(vec[EVN]);
  crash_if_borders_not_allocated(vec[ODD]);
  
  //check buffer size matching
  if(a->buf_size!=nbytes_per_site*bord_vol)
    crash("wrong buffer size (%d) for %d border)",a->buf_size,nbytes_per_site*bord_vol);
  
  //the buffer is lx ordered
  NISSA_PARALLEL_LOOP(ibord_lx,0,bord_vol)
    {
      int dest_lx=loc_vol+ibord_lx;
      int par=loclx_parity[dest_lx];
      int dest_eo=loceo_of_loclx[dest_lx];
      memcpy((char*)(vec[par])+dest_eo*nbytes_per_site,a->recv_buf+ibord_lx*nbytes_per_site,nbytes_per_site);
    }
  
  //we do not sync, because typically we will set borders as valid
}

//start communication using an ev and od border
void spi_start_communicating_ev_and_od_borders(int *nrequest,spi_comm_t *a,void **vec,int nbytes_per_site)
{
  if(!check_borders_valid(vec[EVN])||!check_borders_valid(vec[ODD]))
    {
      GET_THREAD_ID();
      
      //take time and output debugging info
      if(IS_MASTER_THREAD) tot_nissa_comm_time-=take_time();
      verbosity_lv3_master_printf("Starting spi communication of ev or od borders of %s\n",get_vec_name((void*)(*vec)));

      //fill the communicator buffer, start the communication and take time
      fill_spi_sending_buf_with_ev_and_od_vec(a,vec,nbytes_per_site);
      spi_start_comm(a);
      if(IS_MASTER_THREAD) tot_nissa_comm_time+=take_time();
      
      (*nrequest)=8;
    }
  else (*nrequest)=0;
}

//finish communicating
void spi_finish_communicating_ev_and_od_borders(int *nrequest,void **vec,spi_comm_t *a,int nbytes_per_site)
{
  if(((!check_borders_valid(vec[EVN]))||(!check_borders_valid(vec[ODD]))) && (*nrequest==8))
    {
      GET_THREAD_ID();
      
      //take time and make some output
      if(IS_MASTER_THREAD) tot_nissa_comm_time-=take_time();
      verbosity_lv3_master_printf("Finish spi communication of ev or od borders of %s\n",get_vec_name((void*)(*vec)));
	  
      //wait communication to finish, fill back the vector and take time
      spi_comm_wait(a);
      fill_ev_and_od_bord_with_spi_receiving_buf(vec,a,nbytes_per_site);
      if(IS_MASTER_THREAD) tot_nissa_comm_time+=take_time();
      
      set_borders_valid(vec[EVN]);
      set_borders_valid(vec[ODD]);
      
      (*nrequest)=0;
    }
}

//merge the two
void spi_communicate_ev_and_od_borders(void **vec,spi_comm_t *a,int nbytes_per_site)
{
  int nrequest;
  spi_start_communicating_ev_and_od_borders(&nrequest,a,vec,nbytes_per_site);
  spi_finish_communicating_ev_and_od_borders(&nrequest,vec,a,nbytes_per_site);
}
