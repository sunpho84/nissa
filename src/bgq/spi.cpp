#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <malloc.h>
#include <stdlib.h>

#define EXTERN_SPI
#include "spi.hpp"

#include "base/thread_macros.hpp"
#include "base/vectors.hpp"
#include "communicate/communicate.hpp"
#include "routines/ios.hpp"
#include "routines/mpi_routines.hpp"
#ifdef USE_THREADS
 #include "routines/thread.hpp"
#endif

//#define SPI_BARRIER

//#define HINTS_DEBUG

#include <spi/include/kernel/MU.h>
#include <spi/include/mu/InjFifo.h>
#include <spi/include/mu/GIBarrier.h>

namespace nissa
{
  //counter
  volatile uint64_t spi_recv_counter;
  
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
  
  //get the spi coord, grid size and torusicity
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
    
    //get torusicity
    for(int idir=0;idir<5;idir++) spi_dir_is_torus[idir]=ND_GET_TORUS(idir,pers.Network_Config.NetFlags);
    
    //get size
    spi_dir_size[0]=pers.Network_Config.Anodes;
    spi_dir_size[1]=pers.Network_Config.Bnodes;
    spi_dir_size[2]=pers.Network_Config.Cnodes;
    spi_dir_size[3]=pers.Network_Config.Dnodes;
    spi_dir_size[4]=pers.Network_Config.Enodes;
  }
  
  //set the spi neighboirs using MPI communicators
  void set_spi_neighbours()
  {
    //loop on the 8 dirs
    for(int mu=0;mu<4;mu++)
      if(paral_dir[mu])
	for(int bf=0;bf<2;bf++) //backward(0) or forward(1)
	  {
	    //send to one dir, receive from the opposite
	    coords_5D c={0,0,0,0,0};
	    MPI_Sendrecv(spi_rank_coord,sizeof(coords_5D),MPI_BYTE,rank_neigh[bf][mu],0,
			 c,sizeof(coords_5D),MPI_BYTE,rank_neigh[!bf][mu],0,
			 cart_comm,MPI_STATUS_IGNORE);
	    
	    //setup the neighbours and destination
	    MUSPI_SetUpDestination(&spi_neigh[!bf][mu],c[0],c[1],c[2],c[3],c[4]);
	    
	    //copy coords
	    for(int tdir=0;tdir<5;tdir++) spi_dest_coord[(!bf)*4+mu][tdir]=c[tdir];
	  }
  }
  
  //check that the number of hopping to move in each direction is <=1
  void check_all_lattice_neighbours_are_spi_first_neighbours()
  {
    int max_nhop=0;
    for(int mu=0;mu<4;mu++)
      if(paral_dir[mu])
	for(int bf=0;bf<2;bf++)
	  {
	    //compute the number of hop according manhattan metric
	    int nhop=0;
	    for(int idir=0;idir<5;idir++)
	      {
		int coord_neigh=spi_dest_coord[bf*4+mu][idir];
		int off=abs(coord_neigh-spi_rank_coord[idir]);
		if(spi_dir_is_torus[idir]) off=std::min(off,std::min(spi_dir_size[idir]+coord_neigh-spi_rank_coord[idir],
								     spi_dir_size[idir]+spi_rank_coord[idir]-coord_neigh));
		nhop+=off;
	      }
	    max_nhop=std::max(max_nhop,nhop);
	  }
    if(max_nhop>1)
      master_printf("WARNING not all lattice-nodes are first neighbours in SPI grid (non optimal communications)\n");
  }
  
  //find coordinates and set neighbours map
  void set_spi_geometry()
  {
    get_spi_coord();
    set_spi_neighbours();
  }
  
  //set the hints
  void set_spi_hints()
  {
    //database
    const uint64_t fifo_mapP[5]={MUHWI_DESCRIPTOR_TORUS_FIFO_MAP_AP,MUHWI_DESCRIPTOR_TORUS_FIFO_MAP_BP,
				 MUHWI_DESCRIPTOR_TORUS_FIFO_MAP_CP,MUHWI_DESCRIPTOR_TORUS_FIFO_MAP_DP,MUHWI_DESCRIPTOR_TORUS_FIFO_MAP_EP};
    const uint64_t fifo_mapM[5]={MUHWI_DESCRIPTOR_TORUS_FIFO_MAP_AM,MUHWI_DESCRIPTOR_TORUS_FIFO_MAP_BM,
				 MUHWI_DESCRIPTOR_TORUS_FIFO_MAP_CM,MUHWI_DESCRIPTOR_TORUS_FIFO_MAP_DM,MUHWI_DESCRIPTOR_TORUS_FIFO_MAP_EM};
    const uint8_t hintP[4]={MUHWI_PACKET_HINT_AP,MUHWI_PACKET_HINT_BP,MUHWI_PACKET_HINT_CP,MUHWI_PACKET_HINT_DP};
    const uint8_t hintM[4]={MUHWI_PACKET_HINT_AM,MUHWI_PACKET_HINT_BM,MUHWI_PACKET_HINT_CM,MUHWI_PACKET_HINT_DM};
    
    for(int idir=0;idir<8;idir++)
      if(paral_dir[idir%4])
	{
	  //reset fifo map
	  spi_fifo_map[idir]=fifo_mapP[0]|fifo_mapP[1]|fifo_mapP[2]|fifo_mapP[3]|fifo_mapP[4]|
	    fifo_mapM[0]|fifo_mapM[1]|fifo_mapM[2]|fifo_mapM[3]|fifo_mapM[4];
	  spi_hint_ABCD[idir]=spi_hint_E[idir]=0;
	  
#ifdef HINTS_DEBUG
	  printf("----rank %d\n",rank);
#endif
	  //search ABCD hints
	  for(int tdir=0;tdir<4;tdir++)
	    {
	      //check that all but tdir (in ABCD) are equals
	      int ort_eq=1;
	      for(int sdir=0;sdir<5;sdir++)
		if(tdir!=sdir)
		  ort_eq&=(spi_dest_coord[idir][sdir]==spi_rank_coord[sdir]);
	      
#ifdef HINTS_DEBUG
	      char FLAG[]="ABCD";
	      printf("Ext dir rank %d %d, spi dir %d (hint %c), ort_eq %d\n",rank,idir,tdir,FLAG[tdir],ort_eq);
#endif
	      //if are equals, check that we are neighbours
	      if(ort_eq)
		{
		  if((spi_dest_coord[idir][tdir]==spi_rank_coord[tdir]-1)||
		     (spi_dir_is_torus[tdir]&&spi_dest_coord[idir][tdir]==spi_dir_size[tdir]-1&&spi_rank_coord[tdir]==0))
		    {
		      spi_hint_ABCD[idir]=hintM[tdir];
		      spi_fifo_map[idir]=fifo_mapM[tdir];
#ifdef HINTS_DEBUG
		      printf(" rank %d (c %d, n %d, s %d, t %d) found 1\n",rank,spi_rank_coord[tdir],
			     spi_dest_coord[idir][tdir],spi_dir_size[tdir],spi_dir_is_torus[tdir]);
#endif
		    }
		  if((spi_dest_coord[idir][tdir]==spi_rank_coord[tdir]+1)||
		     (spi_dir_is_torus[tdir]&&spi_dest_coord[idir][tdir]==0&&spi_rank_coord[tdir]==spi_dir_size[tdir]-1))
		    {
		      spi_hint_ABCD[idir]=hintP[tdir];
		      spi_fifo_map[idir]=fifo_mapP[tdir];
#ifdef HINTS_DEBUG
		      printf(" rank %d (c %d, n %d, s %d, t %d) found 2\n",rank,spi_rank_coord[tdir],
			     spi_dest_coord[idir][tdir],spi_dir_size[tdir],spi_dir_is_torus[tdir]);
#endif
		    }
		}
	    }
	  
#ifdef HINTS_DEBUG
	  printf("Checking rank %d hint E\n",rank);
#endif
	  //find if all but E are equals
	  int ort_E_eq=1;
	  for(int sdir=0;sdir<4;sdir++) ort_E_eq&=(spi_dest_coord[idir][sdir]==spi_rank_coord[sdir]);
	  
	  //fix E hint
	  if(ort_E_eq)
	    {
	      if(idir<4)
		{
		  spi_hint_E[idir]=MUHWI_PACKET_HINT_EP;
		  spi_fifo_map[idir]=fifo_mapP[4];
#ifdef HINTS_DEBUG
		  printf(" rank %d found 1: %u\n",rank,spi_hint_E[idir]);
#endif
		}
	      else
		{
		  spi_hint_E[idir]=MUHWI_PACKET_HINT_EM;
		  spi_fifo_map[idir]=fifo_mapM[4];
#ifdef HINTS_DEBUG
		  printf(" rank %d found 2: %u\n",rank,spi_hint_E[idir]);
#endif
		}
	    }
	}
    
#ifdef HINTS_DEBUG
    printf("----rank %d\n",rank);
    for(int idir=0;idir<8;idir++)
      printf("Conclusions: rank %d, dir %d, ABCD hint: %d, E: %d\n",rank,idir,spi_hint_ABCD[idir],spi_hint_E[idir]);
    printf("----rank %d\n",rank);
#endif
  }
  
  //initialize the spi communications
  void init_spi()
  {
    //check not to have initialized
    if(!spi_inited)
      {
	verbosity_lv1_master_printf("Starting spi\n");
	
	//check that we do not have more than one process per node
	if(Kernel_ProcessCount()!=1) crash("only one process per node implemented");
	
	//mark as initialized
	spi_inited=true;
	
	//get coordinates, size and rank in the 5D grid
	set_spi_geometry();
	
	//check that all ranks are first neighbours in SPI grid
	check_all_lattice_neighbours_are_spi_first_neighbours();
	
	//allocate bats
	spi_bat_id[0]=0;
	spi_bat_id[1]=1;
	if(Kernel_AllocateBaseAddressTable(0,&spi_bat_gr,2,spi_bat_id,0)) crash("allocating bat");
	
	////////////////////////////////// init the fifos ///////////////////////////////////
	
	//alloc space for the injection fifos
	uint32_t fifo_size=64*NSPI_FIFO;
	for(int ififo=0;ififo<NSPI_FIFO;ififo++) spi_fifo[ififo]=(uint64_t*)memalign(64,fifo_size);
	
	//set default attributes for inj fifo
	Kernel_InjFifoAttributes_t fifo_attrs[NSPI_FIFO];
	memset(fifo_attrs,0,NSPI_FIFO*sizeof(Kernel_InjFifoAttributes_t));
	
	//initialize them with default attributes
	uint32_t fifo_id[NSPI_FIFO];
	for(int ififo=0;ififo<NSPI_FIFO;ififo++) fifo_id[ififo]=ififo;
	if(Kernel_AllocateInjFifos(0,&spi_fifo_sg_ptr,NSPI_FIFO,fifo_id,fifo_attrs)) crash("allocating inj fifos");
	
	//init the MU MMIO for the fifos
	for(int ififo=0;ififo<NSPI_FIFO;ififo++)
	  {
	    //create the memory region
	    Kernel_MemoryRegion_t mem_region;
	    if(Kernel_CreateMemoryRegion(&mem_region,spi_fifo[NSPI_FIFO-1-ififo],fifo_size))
	      crash("creating memory region %d of bytes",ififo,fifo_size);
	    
	    //initialize the fifos
	    if(Kernel_InjFifoInit(&spi_fifo_sg_ptr,fifo_id[ififo],&mem_region,
				  (uint64_t)spi_fifo[NSPI_FIFO-1-ififo]-(uint64_t)mem_region.BaseVa,fifo_size-1))
	      crash("initializing fifo");
	  }
	
	//activate the fifos
	if(Kernel_InjFifoActivate(&spi_fifo_sg_ptr,NSPI_FIFO,fifo_id,KERNEL_INJ_FIFO_ACTIVATE)) crash("activating fifo");
	
	//check alignment
	CRASH_IF_NOT_ALIGNED(recv_buf,64);
	CRASH_IF_NOT_ALIGNED(send_buf,64);
	
	//get physical address of receiving buffer
	Kernel_MemoryRegion_t mem_region;
	if(Kernel_CreateMemoryRegion(&mem_region,recv_buf,recv_buf_size))
	  crash("creating recv_buf memory region of %d bytes",recv_buf_size);
	
	//set the physical address
	if(MUSPI_SetBaseAddress(&spi_bat_gr,spi_bat_id[0],(uint64_t)recv_buf-
				(uint64_t)mem_region.BaseVa+(uint64_t)mem_region.BasePa))
	  crash("setting base address");
	
	//set receive counter bat to MU style atomic PA addr of the receive counter
	if((uint64_t)(&spi_recv_counter)&0x7) crash("recv counter not 8 byte aligned");
	if(Kernel_CreateMemoryRegion(&mem_region,(void*)&spi_recv_counter,sizeof(uint64_t)))
	  crash("creating memory region of %d bytes",sizeof(uint64_t));
	if(MUSPI_SetBaseAddress(&spi_bat_gr,spi_bat_id[1],MUSPI_GetAtomicAddress((uint64_t)&spi_recv_counter-(uint64_t)mem_region.BaseVa+(uint64_t)mem_region.BasePa,MUHWI_ATOMIC_OPCODE_STORE_ADD))) crash("setting base addr");
	
	//reset number of byte to be received
	spi_recv_counter=0;
	
	//get the send buffer physical address
	if(Kernel_CreateMemoryRegion(&mem_region,send_buf,send_buf_size))
	  crash("creating memory region of %d bytes",send_buf_size);
	spi_send_buf_phys_addr=(uint64_t)send_buf-(uint64_t)mem_region.BaseVa+(uint64_t)mem_region.BasePa;
	
	//find hints for descriptors
	set_spi_hints();
	
#ifdef SPI_BARRIER
	//init the barrier
	if(MUSPI_GIBarrierInit(&spi_barrier,0)) crash("initializing the barrier");
#endif
	verbosity_lv2_master_printf("spi initialized\n");
      }
  }
  
  //install spi descriptors, needed for communticaions
  void spi_descriptor_setup(comm_t &in)
  {
    posix_memalign((void**)&in.descriptors,64,8*sizeof(MUHWI_Descriptor_t));
    CRASH_IF_NOT_ALIGNED(in.descriptors,64);
    
    //create one descriptor per direction
    for(int idir=0;idir<8;idir++)
      if(paral_dir[idir%4])
	{
	  //reset the info
	  MUSPI_Pt2PtDirectPutDescriptorInfo_t dinfo;
	  memset(&dinfo,0,sizeof(MUSPI_Pt2PtDirectPutDescriptorInfo_t));
	  
	  //set the parameters
	  dinfo.Base.Payload_Address=spi_send_buf_phys_addr+in.send_offset[idir];
	  dinfo.Base.Message_Length=in.message_length[idir];
	  dinfo.Base.Torus_FIFO_Map=spi_fifo_map[idir];
	  dinfo.Base.Dest=in.spi_dest[idir];
#ifdef HINTS_DEBUG
	  printf("rank %d setting ABCD hint for dir %d to %u\n",rank,idir,spi_hint_ABCD[idir]);
#endif
	  dinfo.Pt2Pt.Hints_ABCD=0; //this is not working...
	  
	  dinfo.Pt2Pt.Misc1=MUHWI_PACKET_USE_DETERMINISTIC_ROUTING|MUHWI_PACKET_DO_NOT_ROUTE_TO_IO_NODE;
	  //|spi_hint_E[idir]; and this as well
#ifdef HINTS_DEBUG
	  printf("rank %d setting E hint for dir %d to %u\n",rank,idir,spi_hint_E[idir]);
#endif
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
      if(paral_dir[idir%4]&&(dir_comm==NULL||dir_comm[idir]))
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
	if(paral_dir[idir%4])
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
  {free(in.descriptors);}
}
