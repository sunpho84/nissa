#ifndef _ALL_TO_ALL_COMM_H
#define _ALL_TO_ALL_COMM_H

#include <map>
#include <utility>
#include <vector>

namespace nissa
{
  struct rank_iel_t
  {
    int rank,iel;
    rank_iel_t(int rank,int iel) : rank(rank),iel(iel) {}
  };
  inline bool operator<(const rank_iel_t &x,const rank_iel_t &y){return x.rank<y.rank||((x.rank==y.rank)&&(x.iel<y.iel));}
  struct gathering_el_t
  {
    rank_iel_t rank_iel_fr;
    int iel_to;
    gathering_el_t(rank_iel_t rank_iel_fr,int iel_to) : rank_iel_fr(rank_iel_fr),iel_to(iel_to) {}
  };
  inline bool operator<(const gathering_el_t &x,const gathering_el_t &y) {return x.rank_iel_fr<y.rank_iel_fr;}
  struct scattering_el_t
  {
    int iel_fr;
    rank_iel_t rank_iel_to;
    scattering_el_t(int iel_fr,rank_iel_t rank_iel_to) : iel_fr(iel_fr),rank_iel_to(rank_iel_to) {}
  };
  struct all_to_all_gathering_list_t : std::vector<gathering_el_t> {};
  struct all_to_all_scattering_list_t : std::vector<scattering_el_t> {};

  struct temp_build_t
  {
    int *nper_rank_to_temp,*nper_rank_fr_temp;
    int *out_buf_cur_per_rank,*in_buf_cur_per_rank;
    std::map<int,int> rank_to_map_list_ranks_to,rank_fr_map_list_ranks_fr;
    temp_build_t()
    {
      nper_rank_to_temp=nissa_malloc("nper_rank_to_temp",nranks,int);
      nper_rank_fr_temp=nissa_malloc("nper_rank_fr_temp",nranks,int);
    }
    ~temp_build_t()
    {
      nissa_free(nper_rank_to_temp);
      nissa_free(nper_rank_fr_temp);
      nissa_free(out_buf_cur_per_rank);
      nissa_free(in_buf_cur_per_rank);
    }
  };

  struct all_to_all_comm_t
  {
    int nel_out,nel_in;
    int nranks_fr,*list_ranks_fr,*in_buf_dest,*nper_rank_fr,*in_buf_off_per_rank;
    int nranks_to,*list_ranks_to,*out_buf_source,*nper_rank_to,*out_buf_off_per_rank;
    
    all_to_all_comm_t(all_to_all_gathering_list_t &gl);
    all_to_all_comm_t(all_to_all_scattering_list_t &sl);
    ~all_to_all_comm_t();
    void communicate(void *out,void *in,int bps,void *buf_out=NULL,void *buf_in=NULL);

    void setup_knowing_where_to_send(all_to_all_scattering_list_t &sl);
    void setup_knowing_what_to_ask(all_to_all_gathering_list_t &gl);
    void setup_nper_rank_ot_temp(int *nper_rank_ot_temp,int *nper_rank_temp);
    void common_setup_part1(temp_build_t &build);
    void common_setup_part2(int nel_note,int *&buf_note,int nranks_note,int *list_ranks_note,int *buf_note_off_per_rank,int *nper_rank_note,int *buf_expl,int nranks_expl,int *list_ranks_expl,int *buf_expl_off_per_rank,int *nper_rank_expl);
    all_to_all_comm_t() {};
  };
}

#endif
