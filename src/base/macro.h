#pragma once

#define EVN 0
#define ODD 1

#define BORDERS_VALID 1
#define EDGES_VALID 2

#define RAD2 1.414213562373095048801688724209l

#define debug_lvl 1

#define nissa_vect_string_length 20
#define nissa_vect_alignment 16

#define nissa_loc_volh_loop(a) for(int a=0;a<loc_volh;a++)
#define nissa_loc_vol_loop(a) for(int a=0;a<loc_vol;a++)

#define ran2_ntab 32                                                                                                         
#define nissa_malloc(a,b,c) nissa_true_malloc(a,b,sizeof(c),#c,__FILE__,__LINE__)
#define nissa_free(a) a=nissa_true_free(a,__FILE__,__LINE__)

#define master_printf(...) master_fprintf(stdout,__VA_ARGS__)
#define debug_fprintf(lvl,fil,...)					\
  if(debug_lvl>=lvl)							\
    {									\
      MPI_Barrier(MPI_COMM_WORLD);					\
      if(rank==0)							\
	{								\
	  fprintf(fil,"Debug message: <<");				\
	  fprintf(fil,__VA_ARGS__);					\
	  fprintf(fil,">> on file %s, line %d\n",__FILE__,__LINE__);	\
	  fflush(fil);							\
	}								\
      MPI_Barrier(MPI_COMM_WORLD);					\
    }

#define debug_printf(lvl,...) debug_fprintf(lvl,stdout,__VA_ARGS__);

#define crash(...) internal_crash(__LINE__,__FILE__,__VA_ARGS__)

#define decript_MPI_error(...) internal_decript_MPI_error(__LINE__,__FILE__,__VA_ARGS__)

#define if_nissa_vect_not_initialized() if(main_nissa_arr!=((char*)&main_nissa_vect)+sizeof(nissa_vect))
