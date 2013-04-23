#ifdef HAVE_CONFIG_H
 #include "config.h"
#endif

#include <mpi.h>

#include "../base/debug.h"
#include "../base/global_variables.h"
#include "../base/vectors.h"
#include "../geometry/geometry_lx.h"
#include "../new_types/new_types_definitions.h"
#include "../routines/ios.h"
#include "../routines/thread.h"

#ifdef BGQ
 #include <malloc.h>
 #include "../bgq/spi.h"
#include <spi/include/mu/Addressing_inlines.h>
#endif

//general set of bufered comm
void install_buffered_comm(buffered_comm_t &in)
{
  //check that buffers are large enough
  if(in.tot_mess_size>nissa_buff_size)
    crash("asking to create a communicator that need %d large buffer (%d allocated)",in.tot_mess_size,nissa_buff_size);
  
  //mark that there is no communication in progress
  in.comm_in_prog=0;
  
  //bgq replacements and original MPI initialization
#ifdef BGQ
  //check alignment
  CRASH_IF_NOT_ALIGNED(nissa_recv_buf,64);
  CRASH_IF_NOT_ALIGNED(nissa_send_buf,64);
  CRASH_IF_NOT_ALIGNED(in.descriptors,64);
  
  //allocate the bat entries
  in.bat_id[0]=spi_nallocated_bat+0;
  in.bat_id[1]=spi_nallocated_bat+1;
  spi_nallocated_bat+=2;
  if(Kernel_AllocateBaseAddressTable(0,&in.spi_bat_gr,2,in.bat_id,0)) crash("allocating bat");
  verbosity_lv3_master_printf("number of spi allocated bat: %d\n",spi_nallocated_bat);
  
  //get physical address of receiving buffer
  Kernel_MemoryRegion_t mem_region;
  if(Kernel_CreateMemoryRegion(&mem_region,nissa_recv_buf,in.tot_mess_size)) crash("creating memory region");
  
  //set the physical address
  if(MUSPI_SetBaseAddress(&in.spi_bat_gr,in.bat_id[0],(uint64_t)nissa_recv_buf-
			  (uint64_t)mem_region.BaseVa+(uint64_t)mem_region.BasePa))
    crash("setting base address");
  
  //set receive counter bat to MU style atomic PA addr of the receive counter
  if((uint64_t)(&in.recv_counter)&0x7) crash("recv counter not 8 byte aligned");
  if(Kernel_CreateMemoryRegion(&mem_region,(void*)&in.recv_counter,sizeof(uint64_t))) crash("creating memory region");  
  if(MUSPI_SetBaseAddress(&in.spi_bat_gr,in.bat_id[1],MUSPI_GetAtomicAddress((uint64_t)&in.recv_counter-(uint64_t)mem_region.BaseVa+(uint64_t)mem_region.BasePa,MUHWI_ATOMIC_OPCODE_STORE_ADD))) crash("setting base addr");
  
  //reset number of byte to be received
  in.recv_counter=0;

  //get the send buffer physical address
  if(Kernel_CreateMemoryRegion(&mem_region,nissa_send_buf,in.tot_mess_size)) crash("creating memory region");
  uint64_t send_buf_phys_addr=(uint64_t)nissa_send_buf-(uint64_t)mem_region.BaseVa+(uint64_t)mem_region.BasePa;
  
  //create one descriptor per direction
  for(int idir=0;idir<8;idir++)
    {
      //reset the info
      MUSPI_Pt2PtDirectPutDescriptorInfo_t dinfo;
      memset(&dinfo,0,sizeof(MUSPI_Pt2PtDirectPutDescriptorInfo_t));
      
      //set the parameters
      dinfo.Base.Payload_Address=send_buf_phys_addr+in.send_offset[idir];
      dinfo.Base.Message_Length=in.message_length[idir];
      dinfo.Base.Torus_FIFO_Map=MUHWI_DESCRIPTOR_TORUS_FIFO_MAP_AM|MUHWI_DESCRIPTOR_TORUS_FIFO_MAP_AP|
	MUHWI_DESCRIPTOR_TORUS_FIFO_MAP_BM|MUHWI_DESCRIPTOR_TORUS_FIFO_MAP_BP|MUHWI_DESCRIPTOR_TORUS_FIFO_MAP_CM|
	MUHWI_DESCRIPTOR_TORUS_FIFO_MAP_CP|MUHWI_DESCRIPTOR_TORUS_FIFO_MAP_DM|MUHWI_DESCRIPTOR_TORUS_FIFO_MAP_DP|
	MUHWI_DESCRIPTOR_TORUS_FIFO_MAP_EM|MUHWI_DESCRIPTOR_TORUS_FIFO_MAP_EP;
      dinfo.Base.Dest=in.spi_dest[idir];
      dinfo.Pt2Pt.Hints_ABCD=0;
      dinfo.Pt2Pt.Misc1=MUHWI_PACKET_USE_DETERMINISTIC_ROUTING|MUHWI_PACKET_DO_NOT_ROUTE_TO_IO_NODE;	
      dinfo.Pt2Pt.Misc2=MUHWI_PACKET_VIRTUAL_CHANNEL_DETERMINISTIC;
      dinfo.Pt2Pt.Skip=8; //for checksumming, skip the header
      dinfo.DirectPut.Rec_Payload_Base_Address_Id=in.bat_id[0];
      dinfo.DirectPut.Rec_Payload_Offset=in.recv_offset[idir];
      dinfo.DirectPut.Rec_Counter_Base_Address_Id=in.bat_id[1];
      dinfo.DirectPut.Rec_Counter_Offset=0;
      dinfo.DirectPut.Pacing=MUHWI_PACKET_DIRECT_PUT_IS_NOT_PACED;
      
      //create the descriptor
      if(MUSPI_CreatePt2PtDirectPutDescriptor(&in.descriptors[idir],&dinfo)) crash("creating the descriptor");
    }
#else
  in.nrequest=0;
  in.imessage=nbuffered_comm_allocated;
#endif
  
  nbuffered_comm_allocated++;
}

//set up a communicator for lx or eo borders
//first 4 communicate to forward nodes, last four to backward nodes
void set_lx_or_eo_buffered_comm(buffered_comm_t &in,int lx_eo,int nbytes_per_site)
{
  int div_coeff=(lx_eo==0)?1:2; //dividing coeff
  
  //copy nbytes_per_site and compute total size
  in.nbytes_per_site=nbytes_per_site;
  in.tot_mess_size=in.nbytes_per_site*bord_vol/div_coeff;
  
  //direction of the halo in receiving node: surface is ordered opposite of halo
  for(int bf=0;bf<2;bf++) 
    for(int mu=0;mu<4;mu++)
      {
	int idir=bf*4+mu;
	
	//set the parameters
	in.send_offset[idir]=(bord_offset[mu]+bord_vol/2*(!bf))*in.nbytes_per_site/div_coeff;
	in.message_length[idir]=bord_dir_vol[mu]*in.nbytes_per_site/div_coeff;
	in.recv_offset[idir]=(bord_offset[mu]+bord_vol/2*bf)*in.nbytes_per_site/div_coeff;
#ifdef BGQ
	in.spi_dest[idir]=spi_neigh[!bf][mu];
#else	
	in.recv_rank[idir]=rank_neigh [bf][mu];
	in.send_rank[idir]=rank_neigh[!bf][mu];
#endif
      }
  
  install_buffered_comm(in);
}

void set_lx_buffered_comm(buffered_comm_t &in,int nbytes_per_site)
{set_lx_or_eo_buffered_comm(in,0,nbytes_per_site);}
void set_eo_buffered_comm(buffered_comm_t &in,int nbytes_per_site)
{set_lx_or_eo_buffered_comm(in,1,nbytes_per_site);}

//wait a communication to finish
void buffered_comm_wait(buffered_comm_t *in)
{
  GET_THREAD_ID();
  
  if(IS_MASTER_THREAD)
    {
      verbosity_lv3_master_printf("Entering spi comm wait\n");
      
      if(in->comm_in_prog)
	{
#ifdef BGQ
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
	  while(in->recv_counter>0)
	    verbosity_lv3_master_printf("%d/%d bytes remaining to be received\n",in->recv_counter,in->tot_mess_size);
#else
	  verbosity_lv3_master_printf("Waiting for %d MPI request\n",in->nrequest);
	  MPI_Waitall(in->nrequest,in->requests,MPI_STATUS_IGNORE);
	  in->comm_in_prog=0;
#endif
	  in->nrequest=0;
	}
      else verbosity_lv3_master_printf("Did not have to wait for any buffered comm\n");
    }
  
  //all threads must wait
  thread_barrier(BUFFERED_COMM_WAIT_BARRIER);
}

//start the communications
void buffered_comm_start(buffered_comm_t *in,int *dir_comm=NULL,int tot_size=-1)
{
  GET_THREAD_ID();
  
  //check that no other buffered communication is in progress
  if(buffered_comm_in_prog) crash("not more than one buffered communication per time possible!");
  
  if(IS_MASTER_THREAD)
    {
#ifdef BGQ
      //reset the counter and wait that all have reset
      spi_global_barrier();
      if(tot_size==-1) in->recv_counter=in->tot_mess_size;
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
#else
      in->nrequest=0;
      
      for(int idir=0;idir<8;idir++)
	if(paral_dir[idir%4] &&(dir_comm==NULL||dir_comm[idir]))
	  {
	    //exchanging the lower surface, from the first half of sending node to the second half of receiving node
	    MPI_Irecv(nissa_recv_buf+in->recv_offset[idir],in->message_length[idir],MPI_CHAR,in->recv_rank[idir],
		      in->imessage,cart_comm,in->requests+(in->nrequest++));
	    MPI_Isend(nissa_send_buf+in->send_offset[idir],in->message_length[idir],MPI_CHAR,in->send_rank[idir],
		      in->imessage,cart_comm,in->requests+(in->nrequest++));
	  }
#endif
      in->comm_in_prog=buffered_comm_in_prog=1;
    }
}

//unset everything (so nothing)
void uninstall_buffered_comm(buffered_comm_t *in)
{
  //wait for any communication to finish
  buffered_comm_wait(in);
}

/////////////////////////////////////// communicating lx vec ///////////////////////////////////

//fill the sending buf using the data inside an lx vec
void fill_buffered_sending_buf_with_lx_vec(buffered_comm_t *a,void *vec)
{
  GET_THREAD_ID();
  
  //check buffer size matching
  if(a->tot_mess_size!=a->nbytes_per_site*bord_vol)
    crash("wrong buffer size (%d) for %d large border)",a->tot_mess_size,a->nbytes_per_site*bord_vol);
  
  //copy one by one the surface of vec inside the sending buffer
  NISSA_PARALLEL_LOOP(ibord,0,bord_vol)
    memcpy(nissa_send_buf+a->nbytes_per_site*ibord,
	   (char*)vec+surflx_of_bordlx[ibord]*a->nbytes_per_site,
	   a->nbytes_per_site);
  
  //wait that all threads filled their portion
  thread_barrier(BUFFERED_COMM_LX_SENDING_BUF_FILL_BARRIER);
}

//extract the information from receiving buffer and put them inside an lx vec
void fill_lx_bord_with_buffered_receiving_buf(void *vec,buffered_comm_t *a)
{
  GET_THREAD_ID();
  
  if(IS_MASTER_THREAD)
    {
      crash_if_borders_not_allocated(vec);
      
      //check buffer size matching
      if(a->tot_mess_size!=a->nbytes_per_site*bord_vol)
	crash("wrong buffer size (%d) for %d large border)",a->tot_mess_size,a->nbytes_per_site*bord_vol);
      
      //the buffer is already ordered as the vec border
      memcpy((char*)vec+loc_vol*a->nbytes_per_site,nissa_recv_buf,a->tot_mess_size);
    }
  
  //we do not sync, because typically we will set borders as valid
}

//start communication using an lx border
void buffered_start_communicating_lx_borders(buffered_comm_t *a,void *vec)
{
  if(!check_borders_valid(vec))
    {
      GET_THREAD_ID();
      
      //take time and write some debug output
      if(IS_MASTER_THREAD) tot_nissa_comm_time-=take_time();
      verbosity_lv3_master_printf("Start buffered communication of lx borders of %s\n",get_vec_name((void*)vec));
      
      //fill the communicator buffer, start the communication and take time
      fill_buffered_sending_buf_with_lx_vec(a,vec);
      buffered_comm_start(a);
      if(IS_MASTER_THREAD) tot_nissa_comm_time+=take_time();
    }
}

//finish communicating
void buffered_finish_communicating_lx_borders(void *vec,buffered_comm_t *a)
{
  if(!check_borders_valid(vec) && (a->nrequest!=0))
    {
      GET_THREAD_ID();
      
      //take note of passed time and write some debug info
      if(IS_MASTER_THREAD) tot_nissa_comm_time-=take_time();
      verbosity_lv3_master_printf("Finish buffered communication of lx borders of %s\n",get_vec_name((void*)vec));

      //wait communication to finish, fill back the vector and take time
      buffered_comm_wait(a);
      fill_lx_bord_with_buffered_receiving_buf(vec,a);
      if(IS_MASTER_THREAD) tot_nissa_comm_time+=take_time();
      
      //set border not valid: this auto sync
      set_borders_valid(vec);
      
      //mark that no communication is in progress
      if(IS_MASTER_THREAD) buffered_comm_in_prog=0;
    }
}

//merge the two
void buffered_communicate_lx_borders(void *vec,buffered_comm_t *a)
{
  if(!check_borders_valid(vec))
    {
      verbosity_lv3_master_printf("Sync buffered communication of lx borders of %s\n",get_vec_name((void*)vec));
      
      buffered_start_communicating_lx_borders(a,vec);
      buffered_finish_communicating_lx_borders(vec,a);
    }
}

/////////////////////////////////////// communicating e/o vec //////////////////////////////////////

//fill the sending buf using the data inside an ev or odd vec
void fill_buffered_sending_buf_with_ev_or_od_vec(buffered_comm_t *a,void *vec,int eo)
{
  GET_THREAD_ID();
  
  //check buffer size matching
  if(a->tot_mess_size!=a->nbytes_per_site*bord_volh)
    crash("wrong buffer size (%d) for %d border)",a->tot_mess_size,a->nbytes_per_site*bord_volh);
  
  //copy one by one the surface of vec inside the sending buffer
  NISSA_PARALLEL_LOOP(ibord,0,bord_volh)
    memcpy(nissa_send_buf+ibord*a->nbytes_per_site,
	   (char*)vec+surfeo_of_bordeo[eo][ibord]*a->nbytes_per_site,a->nbytes_per_site);

  //wait that all threads filled their portion
  thread_barrier(BUFFERED_COMM_EV_OR_OD_SENDING_BUF_FILL_BARRIER);
}

//extract the information from receiving buffer and put them inside an even or odd vec
void fill_ev_or_od_bord_with_buffered_receiving_buf(void *vec,buffered_comm_t *a)
{
  GET_THREAD_ID();
  
  if(IS_MASTER_THREAD)
    {
      crash_if_borders_not_allocated(vec);
      
      //check buffer size matching
      if(a->tot_mess_size!=a->nbytes_per_site*bord_volh)
	crash("wrong buffer size (%d) for %d border)",a->tot_mess_size,a->nbytes_per_site*bord_volh);
      
      //the buffer is already ordered as the vec border
      memcpy((char*)vec+loc_volh*a->nbytes_per_site,nissa_recv_buf,a->tot_mess_size);
    }
  
  //we do not sync, because typically we will set borders as valid
}

//start communication using an ev or od border
void buffered_start_communicating_ev_or_od_borders(buffered_comm_t *a,void *vec,int eo)
{
  if(!check_borders_valid(vec))
    {
      GET_THREAD_ID();
      
      //take time and output debugging info
      if(IS_MASTER_THREAD) tot_nissa_comm_time-=take_time();
      verbosity_lv3_master_printf("Starting buffered communication of ev or od borders of %s\n",get_vec_name((void*)vec));

      //fill the communicator buffer, start the communication and take time
      fill_buffered_sending_buf_with_ev_or_od_vec(a,vec,eo);
      buffered_comm_start(a);
      if(IS_MASTER_THREAD) tot_nissa_comm_time+=take_time();
    }
}

//finish communicating
void buffered_finish_communicating_ev_or_od_borders(void *vec,buffered_comm_t *a)
{
  if(!check_borders_valid(vec) && a->nrequest!=0)
    {
      GET_THREAD_ID();
      
      //take time and make some output
      if(IS_MASTER_THREAD) tot_nissa_comm_time-=take_time();
      verbosity_lv3_master_printf("Finish buffered communication of ev or od borders of %s\n",get_vec_name((void*)vec));
	  
      //wait communication to finish, fill back the vector and take time
      buffered_comm_wait(a);
      fill_ev_or_od_bord_with_buffered_receiving_buf(vec,a);
      if(IS_MASTER_THREAD) tot_nissa_comm_time+=take_time();
      
      set_borders_valid(vec);
    }
}

//merge the two
void buffered_communicate_ev_or_od_borders(void *vec,buffered_comm_t *a,int eo)
{
  if(!check_borders_valid(vec))
    {
      verbosity_lv3_master_printf("Sync buffered communication of ev or od borders of %s\n",get_vec_name((void*)vec));
      
      buffered_start_communicating_ev_or_od_borders(a,vec,eo);
      buffered_finish_communicating_ev_or_od_borders(vec,a);
    }  
}

/////////////////////////////////////// communicating e&o vec //////////////////////////////////////

//fill the sending buf using the data inside an ev and odd vec, using lx style inside buf
void fill_buffered_sending_buf_with_ev_and_od_vec(buffered_comm_t *a,void **vec)
{
  GET_THREAD_ID();
  
  //check buffer size matching
  if(a->tot_mess_size!=a->nbytes_per_site*bord_vol)
    crash("wrong buffer size (%d) for %d border)",a->tot_mess_size,a->nbytes_per_site*bord_vol);
  
  //copy one by one the surface of vec inside the sending buffer
  NISSA_PARALLEL_LOOP(ibord_lx,0,bord_vol)
    {
      //convert lx indexing to eo
      int source_lx=surflx_of_bordlx[ibord_lx];
      int par=loclx_parity[source_lx];
      int source_eo=loceo_of_loclx[source_lx];
      memcpy(nissa_send_buf+ibord_lx*a->nbytes_per_site,(char*)(vec[par])+source_eo*a->nbytes_per_site,
	     a->nbytes_per_site);
    }
  
  //wait that all threads filled their portion
  thread_barrier(BUFFERED_COMM_EV_AND_OD_SENDING_BUF_FILL_BARRIER);
}

//extract the information from receiving buffer and put them inside an even or odd vec
void fill_ev_and_od_bord_with_buffered_receiving_buf(void **vec,buffered_comm_t *a)
{
  GET_THREAD_ID();
  
  //check border allocation
  crash_if_borders_not_allocated(vec[EVN]);
  crash_if_borders_not_allocated(vec[ODD]);
  
  //check buffer size matching
  if(a->tot_mess_size!=a->nbytes_per_site*bord_vol)
    crash("wrong buffer size (%d) for %d border)",a->tot_mess_size,a->nbytes_per_site*bord_vol);
  
  //the buffer is lx ordered
  NISSA_PARALLEL_LOOP(ibord_lx,0,bord_vol)
    {
      int dest_lx=loc_vol+ibord_lx;
      int par=loclx_parity[dest_lx];
      int dest_eo=loceo_of_loclx[dest_lx];
      memcpy((char*)(vec[par])+dest_eo*a->nbytes_per_site,nissa_recv_buf+ibord_lx*a->nbytes_per_site,a->nbytes_per_site);
    }
  
  //we do not sync, because typically we will set borders as valid
}

//start communication using an ev and od border
void buffered_start_communicating_ev_and_od_borders(buffered_comm_t *a,void **vec)
{
  if(!check_borders_valid(vec[EVN])||!check_borders_valid(vec[ODD]))
    {
      GET_THREAD_ID();
      
      //take time and output debugging info
      if(IS_MASTER_THREAD) tot_nissa_comm_time-=take_time();
      verbosity_lv3_master_printf("Starting buffered communication of ev or od borders of %s\n",get_vec_name((void*)(*vec)));

      //fill the communicator buffer, start the communication and take time
      fill_buffered_sending_buf_with_ev_and_od_vec(a,vec);
      buffered_comm_start(a);
      if(IS_MASTER_THREAD) tot_nissa_comm_time+=take_time();
    }
}

//finish communicating
void buffered_finish_communicating_ev_and_od_borders(void **vec,buffered_comm_t *a)
{
  if(a->comm_in_prog && a->nrequest!=0)
    {
      GET_THREAD_ID();
      
      //take time and make some output
      if(IS_MASTER_THREAD) tot_nissa_comm_time-=take_time();
      verbosity_lv3_master_printf("Finish buffered communication of ev or od borders of %s\n",get_vec_name((void*)(*vec)));
	  
      //wait communication to finish, fill back the vector and take time
      buffered_comm_wait(a);
      fill_ev_and_od_bord_with_buffered_receiving_buf(vec,a);
      if(IS_MASTER_THREAD) tot_nissa_comm_time+=take_time();
      
      set_borders_valid(vec[EVN]);
      set_borders_valid(vec[ODD]);
    }
}

//merge the two
void buffered_communicate_ev_and_od_borders(void **vec,buffered_comm_t *a)
{
  buffered_start_communicating_ev_and_od_borders(a,vec);
  buffered_finish_communicating_ev_and_od_borders(vec,a);
}
