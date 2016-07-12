#ifndef _GRID_HPP
#define _GRID_HPP

#include <vector>

#include "geometry/geometry_lx.hpp"

namespace nissa
{
  struct partitioning_t
  {
    std::vector<int> list_fact;
    int icombo,ncombo;
    int factorize_R;
    
    partitioning_t(long long int V,int R);
    
    //! start icombo
    void restart(){icombo=-1;}
    
    //! check not tried all combo
    bool finished(){return icombo>ncombo;}
    
    //! find the partition correseponding to the current combo and check that it satisfies all constraints
    int decrypt_and_validate_partition(coords R_per_dir,coords grid_size,coords min_grid_size,coords fix_R_per_dir);
    
    //! skip all partitioning sharing the passed number of least significative factors
    void skip_combo_of_factors(int ifact);
    
    //find a valid partition
    bool find_next_valid_partition(coords R_per_dir,coords grid_size,coords min_grid_size,coords fix_R_per_dir);
  };
  
  void init_grid(int T,int L);
}

#endif
