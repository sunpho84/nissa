#ifndef _OPTIMIZED_MESONS_2PTS_HPP
#define _OPTIMIZED_MESONS_2PTS_HPP

#include <stdint.h>
#include <stdio.h>
#include <map>
#include <string>
#include <vector>

namespace nissa
{
  ////////////////////// optimized 2 pts //////////////////////
  
  //list of correlation function to which the combo contributes, and weight
  struct corr_id_weight_t: std::map<uint16_t,double> {};
  
  //match of forward and backward components id
  struct forw_back_comp_id_t: std::pair<uint8_t,uint8_t>
    {
      uint8_t &forw() {return this->first;}
      uint8_t &back() {return this->second;}
      forw_back_comp_id_t(uint8_t forw_comp_id,uint8_t back_comp_id)
	{
	  forw()=forw_comp_id;
	  back()=back_comp_id;
	}
    };
  
  //Structure useful to compute a list of correlation function of the type
  // [Re/Im]Tr[G_sink*S_forw*G_sour*G5*S_back^+*G5]
  //we order S as sink-source-re/im.
  //To compute each correlation function we need to summ 32 contribution,
  //relative to the product of the real and imaginary components.
  //Two main patterns are possible, according to the the fact that we must
  //compute the real or imaginary part of the correlation function,
  //and if one of the two dirac matrices are imaginary.
  //We compact the product of the various correlation functions, so that
  //we can compute only once the repeated products.
  struct two_pts_comp_t: std::map<forw_back_comp_id_t,corr_id_weight_t>
    {
      int ncorr;
      void add_sink_source_corr(uint16_t corr_id,double weight,int re_im,uint8_t sink_igamma,uint8_t sour_igamma);
      void print(FILE *fout=stdout);
      void summ_the_loc_forw_back_contractions(double *out,double *S_forw,double *S_back,int nel,int twall);
      void scream();
      std::map<int,std::string> corr_name;
      std::vector<int> pattern_list;
      
      void print_correlations_to_file(FILE *fout,double *corr);
    };
  
  void print_optimized_contractions_to_file(FILE *fout,int ncontr,int *op_sour,int *op_sink,complex *contr,int twall,
					    const char *tag,double norm);
  two_pts_comp_t operator+=(two_pts_comp_t &a,two_pts_comp_t &b);
  two_pts_comp_t operator-=(two_pts_comp_t &a,two_pts_comp_t &b);
  two_pts_comp_t operator*=(two_pts_comp_t &a,double b);
  two_pts_comp_t operator/=(two_pts_comp_t &a,double b);
}

#endif
