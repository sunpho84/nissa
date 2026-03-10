#ifndef _ALL_TO_ALL_COMM_HPP
#define _ALL_TO_ALL_COMM_HPP

#include <map>
#include <vector>

#include "base/vectors.hpp"
#include "geometry/geometry_lx.hpp"

namespace nissa
{
  /// All to all communicators initializing structure
  struct all_to_all_gathering_list_t :
    std::map<int,int>
  {
    int add_conf_link_for_paths(const Coords& g,
				const int& mu);
  };
  
  using all_to_all_scattering_list_t=
    std::vector<std::pair<int,int>>;
  
  struct temp_build_t
  {
    std::vector<int> nper_rank_to_temp;
    
    std::vector<int> nper_rank_fr_temp;
    
    std::vector<int> out_buf_cur_per_rank;
    
    std::vector<int> in_buf_cur_per_rank;
    
    std::map<int,int> rank_to_map_list_ranks_to;
    
    std::map<int,int> rank_fr_map_list_ranks_fr;
    
    temp_build_t();
  };
  
  //all to all communicators
  struct all_to_all_comm_t
  {
    all_to_all_comm_t(const all_to_all_comm_t&)=delete;
    
    int inited{false};
    
    int nel_out{0};
    
    int nel_in{0};
    
    int nranks_fr{0};
    
    std::vector<int> list_ranks_fr;
    
    int* in_buf_dest{nullptr};
    
    std::vector<int> nper_rank_fr;
    
    std::vector<int> in_buf_off_per_rank;
    
    int nranks_to{0};
    
    std::vector<int> list_ranks_to;
    
    int* out_buf_source{nullptr};
    
    std::vector<int> nper_rank_to;
    
    std::vector<int> out_buf_off_per_rank;
    
    all_to_all_comm_t(const all_to_all_gathering_list_t& gl);
    
    all_to_all_comm_t(const all_to_all_scattering_list_t& sl);
    
    all_to_all_comm_t(all_to_all_comm_t&& oth)=default;
    
    ~all_to_all_comm_t()
    {
      destroy();
    }
    
    void destroy()
    {
      if(inited)
	{
	  inited=false;
	  
	  nissa_free(in_buf_dest);
	  nissa_free(out_buf_source);
	}
    }
    
    void communicate(void* out,
		     void* in,
		     const size_t& bps,
		     const int& tag=-1) const;
    
    void setup_knowing_where_to_send(const all_to_all_scattering_list_t& sl);
    
    void setup_knowing_what_to_ask(const all_to_all_gathering_list_t& gl);
    
    void setup_nper_rank_other_temp(std::vector<int>& nper_rank_other_temp,
				    const std::vector<int>& nper_rank_temp);
    
    void common_setup_part1(temp_build_t& build);
    
    void common_setup_part2(const int& nel_note,
			    int*& buf_note,
			    const int& nranks_note,
			    const std::vector<int>& list_ranks_note,
			    const std::vector<int>& buf_note_off_per_rank,
			    const std::vector<int>& nper_rank_note,
			    const std::vector<int>& buf_expl,
			    const int& nranks_expl,
			    const std::vector<int>& list_ranks_expl,
			    const std::vector<int>& buf_expl_off_per_rank,
			    const std::vector<int>& nper_rank_expl);
    
    all_to_all_comm_t() :
      inited(false)
    {
    }
  };
}

#endif
