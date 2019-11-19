#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <mpi.h>
#include <string.h>

#include "base/bench.hpp"
#include "base/debug.hpp"
#include "base/vectors.hpp"
#include "geometry/geometry_eo.hpp"
#include "geometry/geometry_lx.hpp"
#include "geometry/geometry_Leb.hpp"
#include "routines/ios.hpp"
#include "routines/mpi_routines.hpp"
#include "threads/threads.hpp"

#include "communicate.hpp"

#ifdef SPI
 #include <stdlib.h>
 #include "bgq/spi.hpp"
#endif

/*
  general remark: fill the send buf exactly in the same way in which the local border is ordered (bw{0,1,2,3},
  fw{0,1,2,3}, see "communicate.hpp"), because the exchanging routines will automatically take care of reverting
  order, putting up-surface from sending node in dw-border of receiving one. See this 1D example:
  
  send buf N0    send buf N1    ...              recv buf N0    recv buf N1
   ---- ----      ---- ----     ...     --->      ---- ----      ---- ----
  | L0 | H0 |    | L1 | H1 |    ...              | H* | L1 |    | H0 | L2 |
   ---- ----      ---- ----                       ---- ----      ---- ----
*/

namespace nissa
{
  //general set of bufered comm
  void comm_setup(comm_t &comm)
  {
    //mark initialization
    if(comm.initialized) crash("trying to initialize an already initialized communicator!");
    comm.initialized=true;
    
    //check that buffers are large enough
    if(comm.tot_mess_size>std::min(send_buf_size,recv_buf_size))
      crash("asking to create a communicator that need %d large buffer (%d allocated)",
	    comm.tot_mess_size,std::min(send_buf_size,recv_buf_size));
    
    //mark that there is no communication in progress
    comm.comm_in_prog=0;
    
    //bgq replacements and original MPI initialization
#ifdef SPI
    spi_descriptor_setup(comm);
#else
    comm.nrequest=0;
    comm.imessage=ncomm_allocated;
#endif
    
    ncomm_allocated++;
  }
  
  //set up a communicator for lx or eo borders
  //first NDIM communicate to forward nodes, last four to backward nodes
  void set_lx_or_eo_comm(comm_t &comm,int lx_eo,int nbytes_per_site)
  {
    int div_coeff=(lx_eo==0)?1:2; //dividing coeff
    
    //copy nbytes_per_site and compute total size
    comm.nbytes_per_site=nbytes_per_site;
    comm.tot_mess_size=comm.nbytes_per_site*bord_vol/div_coeff;
    
    //direction of the halo in receiving node: surface is ordered opposite of halo
    for(int bf=0;bf<2;bf++)
      for(int mu=0;mu<NDIM;mu++)
	{
	  int idir=bf*NDIM+mu;
	  
	  //set the parameters
	  comm.send_offset[idir]=(bord_offset[mu]+bord_volh*(!bf))*comm.nbytes_per_site/div_coeff;
	  comm.message_length[idir]=bord_dir_vol[mu]*comm.nbytes_per_site/div_coeff;
	  comm.recv_offset[idir]=(bord_offset[mu]+bord_volh*bf)*comm.nbytes_per_site/div_coeff;
#ifndef SPI
	  comm.recv_rank[idir]=rank_neigh [bf][mu];
	  comm.send_rank[idir]=rank_neigh[!bf][mu];
#else
	  comm.spi_dest[idir]=spi_neigh[!bf][mu];
#endif
	}
    
    comm_setup(comm);
  }
  
  void set_lx_comm(comm_t &comm,int nbytes_per_site)
  {set_lx_or_eo_comm(comm,0,nbytes_per_site);}
  void set_eo_comm(comm_t &comm,int nbytes_per_site)
  {set_lx_or_eo_comm(comm,1,nbytes_per_site);}
  
  //check that the communicator is initialized
  void crash_if_not_initialized(comm_t &comm)
  {if(!comm.initialized) crash("using uninitialized communicator!");}
  
  //start the communications
  void comm_start(comm_t &comm,int *dir_comm=NULL,int tot_size=-1)
  {
    GET_THREAD_ID();
    
    //check initialization
    crash_if_not_initialized(comm);
    
    //mark communication as in progress
    comm.comm_in_prog=1;
    
    if(IS_MASTER_THREAD)
      {
#ifdef SPI
	spi_comm_start(comm,dir_comm,tot_size);
#else
	comm.nrequest=0;
	
	for(int idir=0;idir<2*NDIM;idir++)
	  if(paral_dir[idir%NDIM] && (dir_comm==NULL||dir_comm[idir]))
	    {
	      //exchanging the lower surface, from the first half of sending node to the second half of receiving node
	      MPI_Irecv(recv_buf+comm.recv_offset[idir],comm.message_length[idir],MPI_CHAR,comm.recv_rank[idir],
			comm.imessage,cart_comm,comm.requests+(comm.nrequest++));
	      MPI_Isend(send_buf+comm.send_offset[idir],comm.message_length[idir],MPI_CHAR,comm.send_rank[idir],
			comm.imessage,cart_comm,comm.requests+(comm.nrequest++));
	    }
#endif
      }
  }
  
  //wait a communication to finish
  void comm_wait(comm_t &comm)
  {
    GET_THREAD_ID();
    
    //check initialization
    crash_if_not_initialized(comm);
    
    if(IS_MASTER_THREAD)
      {
#ifdef SPI
	verbosity_lv3_master_printf("Entering SPI comm wait\n");
#else
	verbosity_lv3_master_printf("Entering MPI comm wait\n");
#endif
	
	if(comm.comm_in_prog)
	  {
#ifdef SPI
	    spi_comm_wait(comm);
#else
	    verbosity_lv3_master_printf("Waiting for %d MPI request\n",comm.nrequest);
	    MPI_Waitall(comm.nrequest,comm.requests,MPI_STATUS_IGNORE);
#endif
	  }
	else verbosity_lv3_master_printf("Did not have to wait for any buffered comm\n");
      }
    
    //all threads must wait
    THREAD_BARRIER();
    
    //set communications as finished
    comm.comm_in_prog=0;
#ifndef SPI
    comm.nrequest=0;
#endif
  }
  
  //unset everything
  void comm_unset(comm_t &comm)
  {
    //check initialization
    crash_if_not_initialized(comm);
    
    //wait for any communication to finish
    comm_wait(comm);
    
    //mark not initialized
    comm.initialized=false;
    
#ifdef SPI
    spi_descriptor_unset(comm);
#endif
  }
  
  /////////////////////////////////////// communicating lx vec ///////////////////////////////////
  
  //fill the sending buf using the data inside an lx vec
  void fill_sending_buf_with_lx_vec(comm_t &comm,void *vec)
  {
    GET_THREAD_ID();
    
    //check buffer size matching
    if(comm.tot_mess_size!=comm.nbytes_per_site*bord_vol)
      crash("wrong buffer size (%d) for %d large border)",comm.tot_mess_size,comm.nbytes_per_site*bord_vol);
    
    //copy one by one the surface of vec inside the sending buffer
    NISSA_PARALLEL_LOOP(ibord,0,bord_vol)
      memcpy(send_buf+comm.nbytes_per_site*ibord,
	     (char*)vec+surflx_of_bordlx[ibord]*comm.nbytes_per_site,
	     comm.nbytes_per_site);
    NISSA_PARALLEL_LOOP_END;
    
    //wait that all threads filled their portion
    THREAD_BARRIER();
  }
  
  //extract the information from receiving buffer and put them inside an lx vec
  void fill_lx_bord_with_receiving_buf(void *vec,comm_t &comm)
  {
    GET_THREAD_ID();
    
    if(IS_MASTER_THREAD)
      {
	crash_if_borders_not_allocated(vec,comm.nbytes_per_site*(bord_vol+loc_vol));
	
	//check buffer size matching
	if(comm.tot_mess_size!=comm.nbytes_per_site*bord_vol)
	  crash("wrong buffer size (%d) for %d large border)",comm.tot_mess_size,comm.nbytes_per_site*bord_vol);
	
	//the buffer is already ordered as the vec border
	memcpy((char*)vec+loc_vol*comm.nbytes_per_site,recv_buf,comm.tot_mess_size);
      }
    
    //we do not sync, because typically we will set borders as valid
  }
  
  //start communication using an lx border
  void start_communicating_lx_borders(comm_t &comm,void *vec)
  {
    if(!check_borders_valid(vec) && nparal_dir>0)
      {
	GET_THREAD_ID();
	
	//take time and write some debug output
	START_TIMING(tot_comm_time,ntot_comm);
	verbosity_lv3_master_printf("Start communication of lx borders of %s\n",get_vect_name((void*)vec));
	
	//fill the communicator buffer, start the communication and take time
	fill_sending_buf_with_lx_vec(comm,vec);
	comm_start(comm);
	STOP_TIMING(tot_comm_time);
      }
  }
  
  //finish communicating
  void finish_communicating_lx_borders(void *vec,comm_t &comm)
  {
    if(!check_borders_valid(vec) && nparal_dir>0)
      {
	GET_THREAD_ID();
	
	//take note of passed time and write some debug info
	START_TIMING(tot_comm_time,ntot_comm);
	verbosity_lv3_master_printf("Finish communication of lx borders of %s\n",get_vect_name((void*)vec));
	
	//wait communication to finish, fill back the vector and take time
	comm_wait(comm);
	fill_lx_bord_with_receiving_buf(vec,comm);
	STOP_TIMING(tot_comm_time);
	
	//set border not valid: this auto sync
	set_borders_valid(vec);
      }
  }
  
  //merge the two
  void communicate_lx_borders(void *vec,comm_t &comm)
  {
    if(!check_borders_valid(vec))
      {
	verbosity_lv3_master_printf("Sync communication of lx borders of %s\n",get_vect_name((void*)vec));
	
	start_communicating_lx_borders(comm,vec);
	finish_communicating_lx_borders(vec,comm);
      }
  }
  
  /////////////////////////////////////// communicating e/o vec //////////////////////////////////////
  
  //fill the sending buf using the data inside an ev or odd vec
  void fill_sending_buf_with_ev_or_od_vec(comm_t &comm,void *vec,int eo)
  {
    GET_THREAD_ID();
    
    //check buffer size matching
    if(comm.tot_mess_size!=comm.nbytes_per_site*bord_volh)
      crash("wrong buffer size (%d) for %d border)",comm.tot_mess_size,comm.nbytes_per_site*bord_volh);
    
    //copy one by one the surface of vec inside the sending buffer
    NISSA_PARALLEL_LOOP(ibord,0,bord_volh)
      memcpy(send_buf+ibord*comm.nbytes_per_site,
	     (char*)vec+surfeo_of_bordeo[eo][ibord]*comm.nbytes_per_site,comm.nbytes_per_site);
    NISSA_PARALLEL_LOOP_END;
    
    //wait that all threads filled their portion
    THREAD_BARRIER();
  }
  
  //extract the information from receiving buffer and put them inside an even or odd vec
  void fill_ev_or_od_bord_with_receiving_buf(void *vec,comm_t &comm)
  {
    GET_THREAD_ID();
    
    if(IS_MASTER_THREAD)
      {
	crash_if_borders_not_allocated(vec,comm.nbytes_per_site*(bord_volh+loc_volh));
	
	//check buffer size matching
	if(comm.tot_mess_size!=comm.nbytes_per_site*bord_volh)
	  crash("wrong buffer size (%d) for %d border)",comm.tot_mess_size,comm.nbytes_per_site*bord_volh);
	
	//the buffer is already ordered as the vec border
	memcpy((char*)vec+loc_volh*comm.nbytes_per_site,recv_buf,comm.tot_mess_size);
      }
    
    //we do not sync, because typically we will set borders as valid
  }
  
  //start communication using an ev or od border
  void start_communicating_ev_or_od_borders(comm_t &comm,void *vec,int eo)
  {
    if(!check_borders_valid(vec) && nparal_dir>0)
      {
	GET_THREAD_ID();
	
	//take time and output debugging info
	START_TIMING(tot_comm_time,ntot_comm);
	verbosity_lv3_master_printf("Starting communication of ev or od borders of %s\n",get_vect_name((void*)vec));
	
	//fill the communicator buffer, start the communication and take time
	fill_sending_buf_with_ev_or_od_vec(comm,vec,eo);
	comm_start(comm);
	STOP_TIMING(tot_comm_time);
      }
  }
  
  //finish communicating
  void finish_communicating_ev_or_od_borders(void *vec,comm_t &comm)
  {
    if(!check_borders_valid(vec) && nparal_dir>0)
      {
	GET_THREAD_ID();
	
	//take time and make some output
	START_TIMING(tot_comm_time,ntot_comm);
	verbosity_lv3_master_printf("Finish communication of ev or od borders of %s\n",get_vect_name((void*)vec));
	
	//wait communication to finish, fill back the vector and take time
	comm_wait(comm);
	fill_ev_or_od_bord_with_receiving_buf(vec,comm);
	STOP_TIMING(tot_comm_time);
	
	//set border not valid: this auto sync
	set_borders_valid(vec);
      }
  }
  
  //merge the two
  void communicate_ev_or_od_borders(void *vec,comm_t &comm,int eo)
  {
    if(!check_borders_valid(vec) && nparal_dir>0)
      {
	verbosity_lv3_master_printf("Sync communication of ev or od borders of %s\n",get_vect_name((void*)vec));
	
	start_communicating_ev_or_od_borders(comm,vec,eo);
	finish_communicating_ev_or_od_borders(vec,comm);
      }
  }
  
  /////////////////////////////////////// communicating e&o vec //////////////////////////////////////
  
  //fill the sending buf using the data inside an ev and odd vec, using lx style inside buf
  void fill_sending_buf_with_ev_and_od_vec(comm_t &comm,eo_ptr<void> vec)
  {
    GET_THREAD_ID();
    
    //check buffer size matching
    if(comm.tot_mess_size!=comm.nbytes_per_site*bord_vol)
      crash("wrong buffer size (%d) for %d border)",comm.tot_mess_size,comm.nbytes_per_site*bord_vol);
    
    //copy one by one the surface of vec inside the sending buffer
    NISSA_PARALLEL_LOOP(ibord_lx,0,bord_vol)
      {
	//convert lx indexing to eo
	int source_lx=surflx_of_bordlx[ibord_lx];
	int par=loclx_parity[source_lx];
	int source_eo=loceo_of_loclx[source_lx];
	memcpy(send_buf+ibord_lx*comm.nbytes_per_site,(char*)(vec[par])+source_eo*comm.nbytes_per_site,
	       comm.nbytes_per_site);
      }
    NISSA_PARALLEL_LOOP_END;
    
    //wait that all threads filled their portion
    THREAD_BARRIER();
  }
  
  //extract the information from receiving buffer and put them inside an even or odd vec
  void fill_ev_and_od_bord_with_receiving_buf(eo_ptr<void> vec,comm_t &comm)
  {
    GET_THREAD_ID();
    
    //check border allocation
    int min_size=comm.nbytes_per_site*(bord_volh+loc_volh);
    crash_if_borders_not_allocated(vec[EVN],min_size);
    crash_if_borders_not_allocated(vec[ODD],min_size);
    
    //check buffer size matching
    if(comm.tot_mess_size!=comm.nbytes_per_site*bord_vol)
      crash("wrong buffer size (%d) for %d border)",comm.tot_mess_size,comm.nbytes_per_site*bord_vol);
    
    //the buffer is lx ordered
    NISSA_PARALLEL_LOOP(ibord_lx,0,bord_vol)
      {
	int dest_lx=loc_vol+ibord_lx;
	int par=loclx_parity[dest_lx];
	int dest_eo=loceo_of_loclx[dest_lx];
	memcpy((char*)(vec[par])+dest_eo*comm.nbytes_per_site,recv_buf+ibord_lx*comm.nbytes_per_site,comm.nbytes_per_site);
      }
    NISSA_PARALLEL_LOOP_END;
    
    //we do not sync, because typically we will set borders as valid
  }
  
  //start communication using an ev and od border
  void start_communicating_ev_and_od_borders(comm_t &comm,eo_ptr<void> vec)
  {
    if((!check_borders_valid(vec[EVN])||!check_borders_valid(vec[ODD])) && nparal_dir>0)
      {
	GET_THREAD_ID();
	
	//take time and output debugging info
	START_TIMING(tot_comm_time,ntot_comm);
	verbosity_lv3_master_printf("Starting communication of ev and od borders of %s\n",get_vect_name(vec[0]));
	
	//fill the communicator buffer, start the communication and take time
	fill_sending_buf_with_ev_and_od_vec(comm,vec);
	comm_start(comm);
	STOP_TIMING(tot_comm_time);
      }
  }
  
  //finish communicating
  void finish_communicating_ev_and_od_borders(eo_ptr<void> vec,comm_t &comm)
  {
    if(comm.comm_in_prog && nparal_dir>0)
      {
	GET_THREAD_ID();
	
	//take time and make some output
	START_TIMING(tot_comm_time,ntot_comm);
	verbosity_lv3_master_printf("Finish communication of ev and od borders of %s\n",get_vect_name(vec[0]));
	
	//wait communication to finish, fill back the vector and take time
	comm_wait(comm);
	fill_ev_and_od_bord_with_receiving_buf(vec,comm);
	STOP_TIMING(tot_comm_time);
	
	//set border not valid: this auto sync
	set_borders_valid(vec[EVN]);
	set_borders_valid(vec[ODD]);
      }
  }
  
  //merge the two
  void communicate_ev_and_od_borders(eo_ptr<void> vec,comm_t &comm)
  {
    start_communicating_ev_and_od_borders(comm,vec);
    finish_communicating_ev_and_od_borders(vec,comm);
  }
  
  ///////////////////////////////////////////// Lebeo ///////////////////////////////////////////
  
  /////////////////////////////////////// communicating Leblx vec ///////////////////////////////////
  
  //fill the sending buf using the data inside an LebLeblx vec
  void fill_sending_buf_with_Leblx_vec(comm_t &comm,void *vec)
  {
    GET_THREAD_ID();
    
    //check buffer size matching
    if(comm.tot_mess_size!=comm.nbytes_per_site*bord_vol)
      crash("wrong buffer size (%d) for %d large border)",comm.tot_mess_size,comm.nbytes_per_site*bord_vol);
    
    //copy one by one the surface of vec inside the sending buffer
    NISSA_PARALLEL_LOOP(ibord,0,bord_vol)
      memcpy(send_buf+comm.nbytes_per_site*ibord,
	     (char*)vec+surfLeblx_of_bordLeblx[ibord]*comm.nbytes_per_site,
	     comm.nbytes_per_site);
    NISSA_PARALLEL_LOOP_END;
    
    //wait that all threads filled their portion
    THREAD_BARRIER();
  }
  
  //extract the information from receiving buffer and put them inside an Leblx vec
  void fill_Leblx_bord_with_receiving_buf(void *vec,comm_t &comm)
  {fill_lx_bord_with_receiving_buf(vec,comm);}
  
  //start communication using an Leblx border
  void start_communicating_Leblx_borders(comm_t &comm,void *vec)
  {
    if(!check_borders_valid(vec) && nparal_dir>0)
      {
	GET_THREAD_ID();
	
	//take time and write some debug output
	START_TIMING(tot_comm_time,ntot_comm);
	verbosity_lv3_master_printf("Start communication of Leblx borders of %s\n",get_vect_name((void*)vec));
	
	//fill the communicator buffer, start the communication and take time
	fill_sending_buf_with_Leblx_vec(comm,vec);
	comm_start(comm);
	STOP_TIMING(tot_comm_time);
      }
  }
  
  //finish communicating
  void finish_communicating_Leblx_borders(void *vec,comm_t &comm)
  {finish_communicating_lx_borders(vec,comm);}
  
  //merge the two
  void communicate_Leblx_borders(void *vec,comm_t &comm)
  {
    if(!check_borders_valid(vec))
      {
	verbosity_lv3_master_printf("Sync communication of Leblx borders of %s\n",get_vect_name((void*)vec));
	
	start_communicating_Leblx_borders(comm,vec);
	finish_communicating_Leblx_borders(vec,comm);
      }
  }
  
  /////////////////////////////////////// communicating e/o vec //////////////////////////////////////
  
  //fill the sending buf using the data inside an ev or odd vec
  void fill_sending_buf_with_Leb_ev_or_od_vec(comm_t &comm,void *vec,int eo)
  {
    GET_THREAD_ID();
    
    //check buffer size matching
    if(comm.tot_mess_size!=comm.nbytes_per_site*bord_volh)
      crash("wrong buffer size (%d) for %d border)",comm.tot_mess_size,comm.nbytes_per_site*bord_volh);
    
    //copy one by one the surface of vec inside the sending buffer
    NISSA_PARALLEL_LOOP(ibord,0,bord_volh)
      memcpy(send_buf+ibord*comm.nbytes_per_site,
	     (char*)vec+Lebeo_of_loceo[eo][surfeo_of_bordeo[eo][ibord]]*comm.nbytes_per_site,comm.nbytes_per_site);
    NISSA_PARALLEL_LOOP_END;
    
    //wait that all threads filled their portion
    THREAD_BARRIER();
  }
  
  //extract the information from receiving buffer and put them inside an even or odd vec
  void fill_Leb_ev_or_od_bord_with_receiving_buf(void *vec,comm_t &comm)
  {
    GET_THREAD_ID();
    
    if(IS_MASTER_THREAD)
      {
	int min_size=comm.nbytes_per_site*(bord_volh+loc_volh);
	crash_if_borders_not_allocated(vec,min_size);
	
	//check buffer size matching
	if(comm.tot_mess_size!=comm.nbytes_per_site*bord_volh)
	  crash("wrong buffer size (%d) for %d border)",comm.tot_mess_size,comm.nbytes_per_site*bord_volh);
	
	//the buffer is already ordered as the vec border
	memcpy((char*)vec+loc_volh*comm.nbytes_per_site,recv_buf,comm.tot_mess_size);
      }
    
    //we do not sync, because typically we will set borders as valid
  }
  
  //start communication using an ev or od border
  void start_communicating_Leb_ev_or_od_borders(comm_t &comm,void *vec,int eo)
  {
    if(!check_borders_valid(vec) && nparal_dir>0)
      {
	GET_THREAD_ID();
	
	//take time and output debugging info
	START_TIMING(tot_comm_time,ntot_comm);
	verbosity_lv3_master_printf("Starting communication of ev or od borders of %s\n",get_vect_name((void*)vec));
	
	//fill the communicator buffer, start the communication and take time
	fill_sending_buf_with_Leb_ev_or_od_vec(comm,vec,eo);
	comm_start(comm);
	STOP_TIMING(tot_comm_time);
      }
  }
  
  //finish communicating
  void finish_communicating_Leb_ev_or_od_borders(void *vec,comm_t &comm)
  {
    if(!check_borders_valid(vec) && nparal_dir>0)
      {
	GET_THREAD_ID();
	
	//take time and make some output
	START_TIMING(tot_comm_time,ntot_comm);
	verbosity_lv3_master_printf("Finish communication of ev or od borders of %s\n",get_vect_name((void*)vec));
	
	//wait communication to finish, fill back the vector and take time
	comm_wait(comm);
	fill_Leb_ev_or_od_bord_with_receiving_buf(vec,comm);
	STOP_TIMING(tot_comm_time);
	
	//set border not valid: this auto sync
	set_borders_valid(vec);
      }
  }
  
  //merge the two
  void communicate_Leb_ev_or_od_borders(void *vec,comm_t &comm,int eo)
  {
    if(!check_borders_valid(vec) && nparal_dir>0)
      {
	verbosity_lv3_master_printf("Sync communication of ev or od borders of %s\n",get_vect_name((void*)vec));
	
	start_communicating_Leb_ev_or_od_borders(comm,vec,eo);
	finish_communicating_Leb_ev_or_od_borders(vec,comm);
      }
  }
  
  /////////////////////////////////////// communicating e&o vec //////////////////////////////////////
  
  //fill the sending buf using the data inside an ev and odd vec, using Leblx style inside buf
  void fill_sending_buf_with_Leb_ev_and_od_vec(comm_t &comm,eo_ptr<void> vec)
  {
    GET_THREAD_ID();
    
    //check buffer size matching
    if(comm.tot_mess_size!=comm.nbytes_per_site*bord_vol)
      crash("wrong buffer size (%d) for %d border)",comm.tot_mess_size,comm.nbytes_per_site*bord_vol);
    
    //copy one by one the surface of vec inside the sending buffer
    NISSA_PARALLEL_LOOP(ibord_loclx,0,bord_vol)
      {
	//convert Leblx indexing to eo
	int source_loclx=surflx_of_bordlx[ibord_loclx];
	int par=loclx_parity[source_loclx];
	int source_eo=loceo_of_loclx[source_loclx];
	int source_Lebeo=Lebeo_of_loceo[par][source_eo];
	memcpy(send_buf+ibord_loclx*comm.nbytes_per_site,(char*)(vec[par])+source_Lebeo*comm.nbytes_per_site,
	       comm.nbytes_per_site);
      }
    NISSA_PARALLEL_LOOP_END;
    
    //wait that all threads filled their portion
    THREAD_BARRIER();
  }
  
  //extract the information from receiving buffer and put them inside an even or odd vec
  void fill_Leb_ev_and_od_bord_with_receiving_buf(eo_ptr<void> vec,comm_t &comm)
  {fill_ev_and_od_bord_with_receiving_buf(vec,comm);}
  
  //start communication using an ev and od border
  void start_communicating_Leb_ev_and_od_borders(comm_t &comm,eo_ptr<void> vec)
  {
    if((!check_borders_valid(vec[EVN])||!check_borders_valid(vec[ODD])) && nparal_dir>0)
      {
	GET_THREAD_ID();
	
	//take time and output debugging info
	START_TIMING(tot_comm_time,ntot_comm);
	verbosity_lv3_master_printf("Starting communication of Leb ev and od borders of %s\n",get_vect_name((void*)(vec[0])));
	
	//fill the communicator buffer, start the communication and take time
	fill_sending_buf_with_Leb_ev_and_od_vec(comm,vec);
	comm_start(comm);
	STOP_TIMING(tot_comm_time);
      }
  }
  
  //finish communicating
  void finish_communicating_Leb_ev_and_od_borders(eo_ptr<void> vec,comm_t &comm)
  {finish_communicating_ev_and_od_borders(vec,comm);}
  
  //merge the two
  void communicate_Leb_ev_and_od_borders(eo_ptr<void> vec,comm_t &comm)
  {
    start_communicating_Leb_ev_and_od_borders(comm,vec);
    finish_communicating_Leb_ev_and_od_borders(vec,comm);
  }
}
