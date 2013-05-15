#ifdef HAVE_CONFIG_H
 #include "config.h"
#endif

#include <malloc.h>
#include <stdlib.h>

#include <firmware/include/personality.h>
#include <spi/include/mu/Descriptor_inlines.h>
#include <spi/include/mu/Addressing_inlines.h>
#include <spi/include/kernel/MU.h>
#include <spi/include/kernel/location.h>

#include "../src/base/global_variables.h"
#include "../src/base/vectors.h"
#include "../src/routines/ios.h"
#include "../src/routines/thread.h"

//#define SPI_BARRIER

//global barrier for spi
void spi_global_barrier()
{
#ifdef SPI_BARRIER
  if(MUSPI_GIBarrierEnter(&spi_barrier)) crash("while entering spi barrier");
  if(MUSPI_GIBarrierPollWithTimeout(&spi_barrier,60UL*1600000000)) crash("while waiting for barrier to finish");
#else
  MPI_Barrier(MPI_COMM_WORLD);
#endif
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
	coords_5D c={0,0,0,0,0};
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
      
      //allocate bats
      spi_bat_id[0]=0;
      spi_bat_id[1]=1;
      if(Kernel_AllocateBaseAddressTable(0,&spi_bat_gr,2,spi_bat_id,0)) crash("allocating bat");
      verbosity_lv3_master_printf("spi allocated 2 bat\n");
      
      ////////////////////////////////// init the fifos ///////////////////////////////////
      
      //alloc space for the 8 injection fifos
      uint32_t fifo_size=64*8;
      for(int idir=0;idir<8;idir++) spi_fifo[idir]=(uint64_t*)memalign(64,fifo_size);

      //set default attributes for inj fifo
      Kernel_InjFifoAttributes_t fifo_attrs[8];
      for(int idir=0;idir<8;idir++) memset(fifo_attrs,0,8*sizeof(Kernel_InjFifoAttributes_t));
      
      //initialize them with default attributes
      uint32_t fifo_id[8]={0,1,2,3,4,5,6,7};
      if(Kernel_AllocateInjFifos(0,&spi_fifo_sg_ptr,8,fifo_id,fifo_attrs)) crash("allocating inj fifos");
      
      //init the MU MMIO for the fifos
      for(int idir=0;idir<8;idir++)
	{
	  //create the memory region
	  Kernel_MemoryRegion_t mem_region;
	  if(Kernel_CreateMemoryRegion(&mem_region,spi_fifo[7-idir],fifo_size)) crash("creating memory region %d",idir);

	  //initialise the fifos
	  if(Kernel_InjFifoInit(&spi_fifo_sg_ptr,fifo_id[idir],&mem_region,
         	(uint64_t)spi_fifo[7-idir]-(uint64_t)mem_region.BaseVa,fifo_size-1)) crash("initializing fifo");
	}
      
      //activate the fifos
      if(Kernel_InjFifoActivate(&spi_fifo_sg_ptr,8,fifo_id,KERNEL_INJ_FIFO_ACTIVATE)) crash("activating fifo");
  
      //check alignment
      CRASH_IF_NOT_ALIGNED(nissa_recv_buf,64);
      CRASH_IF_NOT_ALIGNED(nissa_send_buf,64);
  
      //get physical address of receiving buffer
      Kernel_MemoryRegion_t mem_region;
      if(Kernel_CreateMemoryRegion(&mem_region,nissa_recv_buf,nissa_buff_size))
	crash("creating nissa_recv_buf memory region");
  
      //set the physical address
      if(MUSPI_SetBaseAddress(&spi_bat_gr,spi_bat_id[0],(uint64_t)nissa_recv_buf-
			      (uint64_t)mem_region.BaseVa+(uint64_t)mem_region.BasePa))
	crash("setting base address");
  
      //set receive counter bat to MU style atomic PA addr of the receive counter
      if((uint64_t)(&spi_recv_counter)&0x7) crash("recv counter not 8 byte aligned");
      if(Kernel_CreateMemoryRegion(&mem_region,(void*)&spi_recv_counter,sizeof(uint64_t)))
	crash("creating memory region");  
      if(MUSPI_SetBaseAddress(&spi_bat_gr,spi_bat_id[1],MUSPI_GetAtomicAddress((uint64_t)&spi_recv_counter-(uint64_t)mem_region.BaseVa+(uint64_t)mem_region.BasePa,MUHWI_ATOMIC_OPCODE_STORE_ADD))) crash("setting base addr");
  
      //reset number of byte to be received
      spi_recv_counter=0;

      //get the send buffer physical address
      if(Kernel_CreateMemoryRegion(&mem_region,nissa_send_buf,nissa_buff_size)) crash("creating memory region");
      spi_send_buf_phys_addr=(uint64_t)nissa_send_buf-(uint64_t)mem_region.BaseVa+(uint64_t)mem_region.BasePa;
      
#ifdef SPI_BARRIER
      //init the barrier
      if(MUSPI_GIBarrierInit(&spi_barrier,0)) crash("initializing the barrier");
#endif
    }
}

/////////////////////////////////////// single communicators ///////////////////////////////////

//install spi descriptors, needed for communticaions
void spi_descriptor_setup(comm_t &in)
{
  posix_memalign((void**)&in.descriptors,64,8*sizeof(MUHWI_Descriptor_t));
  CRASH_IF_NOT_ALIGNED(in.descriptors,64);
  
  //create one descriptor per direction
  for(int idir=0;idir<8;idir++)
    {
      //reset the info
      MUSPI_Pt2PtDirectPutDescriptorInfo_t dinfo;
      memset(&dinfo,0,sizeof(MUSPI_Pt2PtDirectPutDescriptorInfo_t));
      
      //set the parameters
      dinfo.Base.Payload_Address=spi_send_buf_phys_addr+in.send_offset[idir];
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
      dinfo.DirectPut.Rec_Payload_Base_Address_Id=spi_bat_id[0];
      dinfo.DirectPut.Rec_Payload_Offset=in.recv_offset[idir];
      dinfo.DirectPut.Rec_Counter_Base_Address_Id=spi_bat_id[1];
      dinfo.DirectPut.Rec_Counter_Offset=0;
      dinfo.DirectPut.Pacing=MUHWI_PACKET_DIRECT_PUT_IS_NOT_PACED;
      
      //create the descriptor
      if(MUSPI_CreatePt2PtDirectPutDescriptor(&in.descriptors[idir],&dinfo)) crash("creating the descriptor");
    }
}

//start spi communication
void spi_comm_start(comm_t &in,int *dir_comm,int tot_size)
{
  //reset the counter and wait that all have reset
  spi_global_barrier();
  if(tot_size==-1) spi_recv_counter=in.tot_mess_size;
  else spi_recv_counter=tot_size;
  spi_global_barrier();
  
  //start the injection
  for(int idir=0;idir<8;idir++)
    if(dir_comm==NULL||dir_comm[idir])
      {
	verbosity_lv3_master_printf("Injecting %d\n",idir);
	spi_desc_count[idir]=MUSPI_InjFifoInject(MUSPI_IdToInjFifo(idir,&spi_fifo_sg_ptr),&in.descriptors[idir]);
	if(spi_desc_count[idir]>(1ll<<57)) crash("msg_InjFifoInject returned %llu when expecting 1, most likely because there is no room in the fifo",spi_desc_count[idir]);
      }
}

//wait communication end
void spi_comm_wait(comm_t &in)
{
  //wait to send everything
  int wait=1;
  while(wait)
    {
      verbosity_lv3_master_printf("Waiting to finish sending data with spi\n");
      wait=0;
      for(int idir=0;idir<8;idir++)
	wait|=!MUSPI_CheckDescComplete(MUSPI_IdToInjFifo(idir,&spi_fifo_sg_ptr),spi_desc_count[idir]);
    }
  
  spi_global_barrier();
  
  //wait to receive everything
  verbosity_lv3_master_printf("Waiting to finish receiving data with spi\n");
  while(spi_recv_counter>0)
    verbosity_lv3_master_printf("%d/%d bytes remaining to be received\n",spi_recv_counter,in.tot_mess_size);
}

//free communicator (one day unset also them)
void spi_descriptor_unset(comm_t &in)
{
  free(in.descriptors);
}
