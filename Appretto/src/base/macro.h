#pragma once

#define debug_lvl 1

#define appretto_vect_string_length 20

#define appretto_malloc(a,b,c) appretto_true_malloc(a,b,sizeof(c),#c,__FILE__,__LINE__)
#define appretto_free(a) a=appretto_true_free(a,__FILE__,__LINE__)

#define debug(lvl,...) if(debug_lvl>=lvl) {__VA_ARGS__;}
#define master(...) if(rank==0) {__VA_ARGS__;}

#define master_printf(...) master_fprintf(stdout,__VA_ARGS__);
#define debug_fprintf(lvl,fil,...)					\
  {                                                                     \
    MPI_Barrier(MPI_COMM_WORLD);                                        \
    master(debug(lvl,{							\
	  fprintf(fil,"Debug message: <<");				\
	  fprintf(fil,__VA_ARGS__);					\
	  fprintf(fil,">> on file %s, line %d\n",__FILE__,__LINE__);	\
	  fflush(fil);							\
	}))								\
      MPI_Barrier(MPI_COMM_WORLD);					\
  }

#define debug_printf(lvl,...) debug_fprintf(lvl,stdout,__VA_ARGS__);

#define crash(...)                              \
  {                                             \
    debug_fprintf(0,stderr,__VA_ARGS__);	\
    master(MPI_Abort(MPI_COMM_WORLD,1));	\
  }

#define if_appretto_vect_not_initialized() if(main_appretto_arr!=((char*)&main_appretto_vect)+sizeof(appretto_vect))
