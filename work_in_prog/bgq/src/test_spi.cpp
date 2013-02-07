#include "nissa.h"

#include <mpi.h>

//#include <hwi/include/bqc/A2_core.h>
//#include <hwi/include/bqc/A2_inlines.h>
//#include <hwi/include/bqc/MU_PacketCommon.h>
#include <firmware/include/personality.h>
//#include <spi/include/mu/Descriptor.h>
#include <spi/include/mu/Descriptor_inlines.h>
//#include <spi/include/mu/InjFifo.h>
//#include <spi/include/mu/Addressing.h>
//#include <spi/include/mu/Addressing_inlines.h>
//#include <spi/include/mu/GIBarrier.h>
#include <spi/include/kernel/MU.h>
//#include <spi/include/kernel/process.h>
#include <spi/include/kernel/location.h>

//the eights directions are: t-,t+,x-,x+,etc

//type to hold the 5D coordinates
typedef uint8_t coords_5D[5];

//spi rank, coordinates and size of 5D grid
int spi_rank;
coords_5D spi_rank_coord,nspi_rank_dir;

//neighbours int the 4 dirs
MUHWI_Destination_t spi_neigh[2][4];

//compute the spi rank of passed coord
int spi_rank_of_coords(coords_5D c)
{
  //multiply shift
  int out=0;
  for(int i=0;i<5;i++)
    out=out*nspi_rank_dir[i]+c[i];

  return out;
}

//set the spi rank
void set_spi_coord_size_rank(Personality_t &pers)
{
  //copy the coords
  spi_rank_coord[0]=pers.Network_Config.Acoord;
  spi_rank_coord[1]=pers.Network_Config.Bcoord;
  spi_rank_coord[2]=pers.Network_Config.Ccoord;
  spi_rank_coord[3]=pers.Network_Config.Dcoord;
  spi_rank_coord[4]=pers.Network_Config.Ecoord;

  //copy the size
  nspi_rank_dir[0]=pers.Network_Config.Anodes;
  nspi_rank_dir[1]=pers.Network_Config.Bnodes;
  nspi_rank_dir[2]=pers.Network_Config.Cnodes;
  nspi_rank_dir[3]=pers.Network_Config.Dnodes;
  nspi_rank_dir[4]=pers.Network_Config.Enodes;

  //compute the local spi rank
  spi_rank=spi_rank_of_coords(spi_rank_coord);
  
  //print out the spi grid size
  verbosity_lv2_master_printf("Spi grid: %d ",nspi_rank_dir[0]);
  for(int i=1;i<5;i++)
    verbosity_lv2_master_printf("x %d",nspi_rank_dir[i]);
  verbosity_lv2_master_printf("\n");
}

//set the spi neighboirs using MPI communicators
void set_spi_neighbours()
{
  //loop on the 8 dirs
  for(int mu=0;mu<4;mu++)
    for(int bf=0;bf<2;bf++) //backward(0) or forward(1)
      {
	//send to one dir, receive from the other
	coords_5D c;
	MPI_Sendrecv(spi_rank_coord,sizeof(coords_5D),MPI_BYTE,rank_neigh[bf][mu],0,
		     c,sizeof(coords_5D),MPI_BYTE,rank_neigh[!bf][mu],0,
		     cart_comm,MPI_STATUS_IGNORE);
	//setup what just received
	MUSPI_SetUpDestination(&spi_neigh[!bf][mu],c[0],c[1],c[2],c[3],c[4]);
      }
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
      
      //get the personality
      Personality_t pers;
      Kernel_GetPersonality(&pers,sizeof(pers));
      
      //get coordinates, size and rank in the 5D grid
      set_spi_coord_size_rank(pers);
      
      //set the spi neighbours
      set_spi_neighbours();
    }
}

int main(int narg,char **arg)
{
  init_nissa(narg,arg);
  
  init_grid(16,16);
  
  init_spi();
  
  close_nissa();
  
  return 0;
}
