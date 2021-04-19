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
  /// Total number of ranks
  ///
  /// Internal implementation
  EXTERN_RANK int _nranks;
  
  /// Total number of ranks
  inline constexpr const int& nranks=_nranks;
  
  /// This rank
  ///
  /// Internal implementation
  EXTERN_RANK int _rank;
  
  /// This rank
  inline constexpr const int& rank=_rank;
  
  /// Master rank which prints
  ///
  /// Internal implementation
  EXTERN_RANK int _master_rank INIT_RANK_TO(=0);
  
  /// Master rank which prints
  inline constexpr const int& master_rank=_master_rank;
  
  /// Cartesian rank. Is this used anywhere?
  ///
  /// Internal implementation
  EXTERN_RANK int _cart_rank;
  
  /// Cartesian rank. Is this used anywhere?
  inline constexpr const int& cart_rank=_cart_rank;
  
  /// Return whether this is master rank
  inline bool is_master_rank()
  {
    return rank==master_rank;
  }
}

#undef EXTERN_RANK
#undef INIT_RANK_TO

#endif
