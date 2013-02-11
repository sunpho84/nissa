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

//get the spi coord and grid size
void get_SPI_coord()
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
void set_SPI_neighbours()
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

//initialize the spi communications
void init_SPI()
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
      get_SPI_coord();
      
      //set the spi neighbours
      set_SPI_neighbours();
      
      ////////////////////////////////// init the fifos ///////////////////////////////////
      
      //alloc space for the 8 injection fifos
      uint32_t fifo_size=64*8;
      for(int idir=0;idir<8;idir++)
	spi_fifo[idir]=(uint64_t*)memalign(64,fifo_size);
      
      //set default attributes for inj fifo
      Kernel_InjFifoAttributes_t fifo_attrs[8];
      for(int idir=0;idir<8;idir++)
	memset(&fifo_attrs[idir],0,sizeof(Kernel_InjFifoAttributes_t));
      
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
}

struct spi_buff_t
{
  uint64_t buf_size;
  char *recv_buf;
  char *send_buf;
  volatile uint64_t recv_counter;
  uint64_t send_buf_phys_addr;
};

//create the base address table
void create_SPI_bat(spi_buff_t &in)
{
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
  
  //get the send buffer physical address
  if(Kernel_CreateMemoryRegion(&mem_region,in.send_buf,in.buf_size)) crash("creating memory region");
  in.send_buf_phys_addr=(uint64_t)in.send_buf-(uint64_t)mem_region.BaseVa+(uint64_t)mem_region.BasePa;
}

//allocate the buffers
void create_SPI_buffers(spi_buff_t &in,int buf_size)
{
  //copy sizeand allocate
  in.buf_size=buf_size;
  in.recv_buf=(char*)memalign(64,buf_size);
  in.send_buf=(char*)memalign(64,buf_size);
}

void test_SPI_comm()
{
  spi_buff_t a;
  create_SPI_buffers(a,64);
  create_SPI_bat(a);
  master_printf("test passed\n");
}

//create the descriptors
//spi_descriptors=(char*)memalign(64,sizeof(MUHWI_Descriptor_t)*8);
  
