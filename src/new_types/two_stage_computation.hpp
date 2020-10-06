#ifndef _TWO_STAGE_COMPUTATION_HPP
#define _TWO_STAGE_COMPUTATION_HPP

#include "base/vectors.hpp"

namespace nissa
{
  // Type to hold the position of linear algebra output data (see "two_stage_computations" doc for explenations)
  struct two_stage_computation_pos_t
  {
    int *inter_fr_in_pos; //offset for intermediate result
    int *final_fr_inter_pos; //offset for final result from intermediate
    int *inter_fr_recv_pos; //offset for intermediate from nissa_recv_buf
    two_stage_computation_pos_t()
    {
      inter_fr_in_pos=final_fr_inter_pos=inter_fr_recv_pos=NULL;
    }
    void free()
    {
      if(inter_fr_in_pos!=NULL) nissa_free(inter_fr_in_pos);
      if(final_fr_inter_pos!=NULL) nissa_free(final_fr_inter_pos);
      if(inter_fr_recv_pos!=NULL) nissa_free(inter_fr_recv_pos);
    }
    ~two_stage_computation_pos_t() {free();}
  };
}

#endif
