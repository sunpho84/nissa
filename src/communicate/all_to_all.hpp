#ifndef _ALL_TO_ALL_COMM_HPP
#define _ALL_TO_ALL_COMM_HPP

#include <algorithm>
#include <map>
#include <vector>

#include "base/vectors.hpp"
#include "geometry/geometry_lx.hpp"

namespace nissa
{
  //all to all communicators initializing structure
  struct all_to_all_gathering_list_t : std::map<int,int>
  {
    int add_conf_link_for_paths(const coords_t& g,const int& mu);
  };
  
  using all_to_all_scattering_list_t=std::vector<std::pair<int,int>>;
  
  struct temp_build_t
  {
    int *nper_rank_to_temp,*nper_rank_fr_temp;
    int *out_buf_cur_per_rank,*in_buf_cur_per_rank;
    std::map<int,int> rank_to_map_list_ranks_to,rank_fr_map_list_ranks_fr;
    temp_build_t();
    ~temp_build_t();
  };
  
  //all to all communicators
  struct all_to_all_comm_t
  {
  private:
    all_to_all_comm_t(const all_to_all_comm_t&)=delete;
  public:
    int inited{false};
    int nel_out{0},nel_in{0};
    int nranks_fr{0},*list_ranks_fr{nullptr},*in_buf_dest{nullptr},*nper_rank_fr{nullptr},*in_buf_off_per_rank{nullptr};
    int nranks_to{0},*list_ranks_to{nullptr},*out_buf_source{nullptr},*nper_rank_to{nullptr},*out_buf_off_per_rank{nullptr};
    
    all_to_all_comm_t(const all_to_all_gathering_list_t &gl);
    all_to_all_comm_t(const all_to_all_scattering_list_t &sl);
    all_to_all_comm_t(all_to_all_comm_t &&oth)
    {
      std::swap(inited,oth.inited);
      std::swap(nel_out,oth.nel_out);
      std::swap(nel_in,oth.nel_in);
      std::swap(nranks_fr,oth.nranks_fr);
      std::swap(list_ranks_fr,oth.list_ranks_fr);
      std::swap(in_buf_dest,oth.in_buf_dest);
      std::swap(nper_rank_fr,oth.nper_rank_fr);
      std::swap(in_buf_off_per_rank,oth.in_buf_off_per_rank);
      std::swap(nranks_to,oth.nranks_to);
      std::swap(list_ranks_to,oth.list_ranks_to);
      std::swap(out_buf_source,oth.out_buf_source);
      std::swap(nper_rank_to,oth.nper_rank_to);
      std::swap(out_buf_off_per_rank,oth.out_buf_off_per_rank);
    }
    
    ~all_to_all_comm_t(){destroy();}
    void destroy()
    {
      if(inited)
	{
	  inited=false;
	  
	  nissa_free(list_ranks_to);
	  nissa_free(list_ranks_fr);
	  nissa_free(in_buf_dest);
	  nissa_free(out_buf_source);
	  nissa_free(nper_rank_fr);
	  nissa_free(nper_rank_to);
	  nissa_free(out_buf_off_per_rank);
	  nissa_free(in_buf_off_per_rank);
	}
    }
    
    void communicate(void *out,
		     void *in,
		     size_t bps,
		     void *buf_out=NULL,
		     void *buf_in=NULL,
		     int tag=-1) const;
    
    void setup_knowing_where_to_send(const all_to_all_scattering_list_t &sl);
    void setup_knowing_what_to_ask(const all_to_all_gathering_list_t &gl);
    void setup_nper_rank_other_temp(int *nper_rank_other_temp,int *nper_rank_temp);
    void common_setup_part1(temp_build_t &build);
    void common_setup_part2(int nel_note,int *&buf_note,int nranks_note,int *list_ranks_note,int *buf_note_off_per_rank,int *nper_rank_note,int *buf_expl,int nranks_expl,int *list_ranks_expl,int *buf_expl_off_per_rank,int *nper_rank_expl);
    all_to_all_comm_t() : inited(false) {};
  };
}

#endif
