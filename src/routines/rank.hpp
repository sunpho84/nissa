#ifndef _RANK_HPP
#define _RANK_HPP

#ifndef EXTERN_RANK
# define EXTERN_RANK extern
# define INIT_RANK_TO(...)
#else
# define INIT_RANK_TO(ARGS...) ARGS
#endif

namespace nissa
{
  //ranks
  EXTERN_RANK int rank,nranks,cart_rank;
  
  EXTERN_RANK int master_rank INIT_RANK_TO(=0);
  
  /// Return whether this is master rank
  inline bool is_master_rank()
  {
    return rank==master_rank;
  }
}

#undef EXTERN_RANK
#undef INIT_RANK_TO

#endif
