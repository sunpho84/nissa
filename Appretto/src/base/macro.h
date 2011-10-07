#pragma once

#define RAD2 1.414213562373095048801688724209l

#define debug_lvl 1

#define appretto_vect_string_length 20

#define ran2_ntab 32                                                                                                         
#define appretto_malloc(a,b,c) appretto_true_malloc(a,b,sizeof(c),#c,__FILE__,__LINE__)
#define appretto_free(a) a=appretto_true_free(a,__FILE__,__LINE__)

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

#define if_appretto_vect_not_initialized() if(main_appretto_arr!=((char*)&main_appretto_vect)+sizeof(appretto_vect))
