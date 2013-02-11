#include "nissa.h"

#include "add_var.h"

#include <malloc.h>

//#include <hwi/include/bqc/A2_core.h>
//#include <hwi/include/bqc/A2_inlines.h>
//#include <hwi/include/bqc/MU_PacketCommon.h>
#include <firmware/include/personality.h>
//#include <spi/include/mu/Descriptor.h>
#include <spi/include/mu/Descriptor_inlines.h>
//#include <spi/include/mu/InjFifo.h>
//#include <spi/include/mu/Addressing.h>
#include <spi/include/mu/Addressing_inlines.h>
//#include <spi/include/mu/GIBarrier.h>
#include <spi/include/kernel/MU.h>
//#include <spi/include/kernel/process.h>
#include <spi/include/kernel/location.h>

//global barrier for spi
void spi_global_barrier()
{
  if(MUSPI_GIBarrierEnter(&spi_barrier)) crash("while entering spi barrier");
  if(MUSPI_GIBarrierPollWithTimeout(&spi_barrier,60UL*1600000000))
    {
      DelayTimeBase(200000000000UL);
      crash("while waiting for barrier to finish");
    }
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
  //flag for spi initialization
  static int nissa_spi_inited=false;
  
  //check not to have initialized
  if(!nissa_spi_inited)
    {
      //check that we do not have more than one process per node
      if(Kernel_ProcessCount()!=1) crash("Only one process per node implemented!");
      
      //mark as initialized
      nissa_spi_inited=true;
      
      //get coordinates, size and rank in the 5D grid
      set_spi_geometry();
      
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

	  //initialise the Fifos
	  if(Kernel_InjFifoInit(&spi_fifo_sg_ptr,fifo_id[idir],&mem_region,
				(uint64_t)spi_fifo[7-idir]-(uint64_t)mem_region.BaseVa,fifo_size-1)) crash("initializing fifo");
	}
      
      //activate the fifos
      if(Kernel_InjFifoActivate(&spi_fifo_sg_ptr,8,fifo_id,KERNEL_INJ_FIFO_ACTIVATE)) crash("activating fifo");
    }
  
  //init the barrier
  if(MUSPI_GIBarrierInit(&spi_barrier,0)) crash("initialing the barrier");
}

//create the spi buffer
void allocate_spi_comm(spi_comm_t &in,int buf_size)
{
  //reset comm in progress
  in.comm_in_prog=0;
  
  /////////////////////////////////////// allocate buffers ///////////////////////////////////////
  
  in.buf_size=buf_size;
  in.recv_buf=(char*)memalign(64,buf_size);
  in.send_buf=(char*)memalign(64,buf_size);
  
  //////////////////////////////// allocate base address table (bat) /////////////////////////////
  
  //allocate the bat entries
  uint32_t bat_id[2]={0,1};
  MUSPI_BaseAddressTableSubGroup_t spi_bat_gr;
  if(Kernel_AllocateBaseAddressTable(0,&spi_bat_gr,2,bat_id,0)) crash("allocating bat");
				   
  //get physical address of receiving buffer
  Kernel_MemoryRegion_t mem_region;
  if(Kernel_CreateMemoryRegion(&mem_region,in.recv_buf,in.buf_size)) crash("creating memory region");
  
  //set the physical address
  if(MUSPI_SetBaseAddress(&spi_bat_gr,0,(uint64_t)in.recv_buf-(uint64_t)mem_region.BaseVa+(uint64_t)mem_region.BasePa))
    crash("setting base address");
  
  //set receive counter bat to MU style atomic PA addr of the receive counter
  if((uint64_t)(&in.recv_counter)&0x7) crash("recv counter not 8 byte aligned");
  if(Kernel_CreateMemoryRegion(&mem_region,(void*)&in.recv_counter,sizeof(uint64_t))) crash("creating memory region");
  if(MUSPI_SetBaseAddress(&spi_bat_gr,1,
			  MUSPI_GetAtomicAddress((uint64_t)&in.recv_counter-(uint64_t)mem_region.BaseVa+(uint64_t)mem_region.BasePa,
						 MUHWI_ATOMIC_OPCODE_STORE_ADD))) crash("setting base addr");
}

//set up a communicator for lx borders
void set_lx_spi_comm(spi_comm_t &in,int bytes_per_site)
{
  //allocate in and out buffers
  allocate_spi_comm(in,bytes_per_site*bord_vol);
  
  //get the send buffer physical address
  Kernel_MemoryRegion_t mem_region;
  if(Kernel_CreateMemoryRegion(&mem_region,in.send_buf,in.buf_size)) crash("creating memory region");
  uint64_t send_buf_phys_addr=(uint64_t)in.send_buf-(uint64_t)mem_region.BaseVa+(uint64_t)mem_region.BasePa;
  
  //create one descriptor per direction
  int idir=0;
  uint64_t offset=0;
  for(int bf=0;bf<2;bf++) //direction of the halo in receiving node: surface is ordered as halo
    for(int mu=0;mu<4;mu++)
      {
	//reset the info
	MUSPI_Pt2PtDirectPutDescriptorInfo_t dinfo;
	memset(&dinfo,0,sizeof(MUSPI_Pt2PtDirectPutDescriptorInfo_t));
	
	//set the parameters
	dinfo.Base.Payload_Address=send_buf_phys_addr+offset;
	dinfo.Base.Message_Length=bord_dir_vol[mu]*bytes_per_site;
	dinfo.Base.Torus_FIFO_Map=MUHWI_DESCRIPTOR_TORUS_FIFO_MAP_AM|MUHWI_DESCRIPTOR_TORUS_FIFO_MAP_AP|
	  MUHWI_DESCRIPTOR_TORUS_FIFO_MAP_BM|MUHWI_DESCRIPTOR_TORUS_FIFO_MAP_BP|MUHWI_DESCRIPTOR_TORUS_FIFO_MAP_CM|
	  MUHWI_DESCRIPTOR_TORUS_FIFO_MAP_CP|MUHWI_DESCRIPTOR_TORUS_FIFO_MAP_DM|MUHWI_DESCRIPTOR_TORUS_FIFO_MAP_DP|
	  MUHWI_DESCRIPTOR_TORUS_FIFO_MAP_EM|MUHWI_DESCRIPTOR_TORUS_FIFO_MAP_EP;
	dinfo.Base.Dest=spi_neigh[!bf][mu];
	dinfo.Pt2Pt.Misc1=MUHWI_PACKET_USE_DETERMINISTIC_ROUTING|MUHWI_PACKET_DO_NOT_ROUTE_TO_IO_NODE;	
	dinfo.Pt2Pt.Misc2=MUHWI_PACKET_VIRTUAL_CHANNEL_DETERMINISTIC;
	dinfo.Pt2Pt.Skip=8; //for checksumming, skip the header
	dinfo.DirectPut.Rec_Payload_Base_Address_Id=0;
	dinfo.DirectPut.Rec_Payload_Offset=offset;
	dinfo.DirectPut.Rec_Counter_Base_Address_Id=1;
	dinfo.DirectPut.Rec_Counter_Offset=0;
	dinfo.DirectPut.Pacing=MUHWI_PACKET_DIRECT_PUT_IS_NOT_PACED;
	
	//increment the offset
	offset+=bord_dir_vol[mu]*bytes_per_site;

	//create the descriptor
	if(MUSPI_CreatePt2PtDirectPutDescriptor(&in.descriptors[idir],&dinfo)) crash("creating the descriptor");
	
	idir++;
      }
}

//wait a communication to finish
void spi_comm_wait(spi_comm_t &in)
{
  //wait to recive everything
  if(in.comm_in_prog)
    while(in.recv_counter>0)
      master_printf("%llu bytes still to be received\n",in.recv_counter);
  
  //wait to send everything
  while(in.comm_in_prog)
    for(int idir=0;idir<8;idir++)
      {
	int o=MUSPI_CheckDescComplete(MUSPI_IdToInjFifo(idir,&spi_fifo_sg_ptr),spi_desc_count[idir]);
	in.comm_in_prog|=!o;
	master_printf("o[%d]=%d\n",idir,o);
      }
}

//start sending
void spi_start_sending(spi_comm_t &in)
{
  //wait for any previous communication to finish and mark as new started 
  spi_comm_wait(in);
  in.comm_in_prog=1;
  
  //reset the counter
  in.recv_counter=in.buf_size;
  
  //wait that all have resetted
  spi_global_barrier();
  
  //start the injection
  for(int idir=0;idir<8;idir++)
    {
      spi_desc_count[idir]=MUSPI_InjFifoInject(MUSPI_IdToInjFifo(idir,&spi_fifo_sg_ptr),in.descriptors);
      if(spi_desc_count[idir]==-1) crash("msg_InjFifoInject failed, most likely because there is no room in the fifo");
    }
}

void test_spi_comm()
{
  spi_comm_t a;
  set_lx_spi_comm(a,sizeof(double));
  
  //fill in the surface the rank
  for(int ivol=0;ivol<bord_vol;ivol++)
    ((double*)a.send_buf)[ivol]=rank;
  
  spi_start_sending(a);
  spi_comm_wait(a);
  
  master_printf("test passed\n");
}
