#include <nissa.hpp>

#define EXTERN_PARS
#include "pars.hpp"

#include "prop.hpp"
#include "conf.hpp"
#include "contr.hpp"

namespace nissa
{
  //read all photon pars
  void read_photon_pars()
  {
    //Zero mode subtraction
    char zero_mode_sub_str[100];
    read_str_str("ZeroModeSubtraction",zero_mode_sub_str,100);
    
    if(strncasecmp(zero_mode_sub_str,"PECIONA",100)==0) photon.zms=PECIONA;
    else
      if(strncasecmp(zero_mode_sub_str,"UNNO_ALEMANNA",100)==0) photon.zms=UNNO_ALEMANNA;
      else
	if(strncasecmp(zero_mode_sub_str,"ONLY_100",100)==0) photon.zms=ONLY_100;
	else crash("Unkwnown zero mode subtraction: %s",zero_mode_sub_str);
    
    //gauge for photon propagator
    char photon_gauge_str[100];
    read_str_str("PhotonGauge",photon_gauge_str,100);
    if(strncasecmp(photon_gauge_str,"FEYNMAN",100)==0) photon.alpha=FEYNMAN_ALPHA;
    else
      if(strncasecmp(photon_gauge_str,"LANDAU",100)==0) photon.alpha=LANDAU_ALPHA;
      else
	if(strncasecmp(photon_gauge_str,"LANDAU",100)==0) read_str_double("Alpha",&photon.alpha);
	else crash("Unkwnown photon gauge: %s",photon_gauge_str);
    
    //discretization for photon propagator
    char photon_discrete_str[100];
    read_str_str("PhotonDiscretization",photon_discrete_str,100);
    if(strncasecmp(photon_discrete_str,"WILSON",100)==0) photon.c1=WILSON_C1;
    else
      if(strncasecmp(photon_discrete_str,"TLSYM",100)==0) photon.c1=TLSYM_C1;
      else crash("Unkwnown photon discretization: %s",photon_discrete_str);
    
    //compute the tadpole summing all momentum
    compute_tadpole(tadpole,photon);
  }
  
  //meson tags
  const int nmes2pts_known=20;
  enum mes2pts_known_t                      { P5P5 , GIP5 , P5GI , V0V0 , AKAK , VKVK , VKTK , TKVK , TKTK , BKBK , GIS0 , S0GI , V0P5 , VKP5 , S0S0 , A0A0 , AKBK , BKAK , V0S0 , S0V0};
  const char mes2pts_tag[nmes2pts_known][5]={"P5P5","GIP5","P5GI","V0V0","AKAK","VKVK","VKTK","TKVK","TKTK","BKBK","GIS0","S0GI","V0P5","VKP5","S0S0","A0A0","AKBK","BKAK","V0S0","S0V0"};
  mes2pts_known_t read_2pts_tag()
  {
    //read the tag
    char tag[10];
    read_str(tag,10);
    
    //convert to int
    int out=0;
    while(strcasecmp(tag,mes2pts_tag[out]) and out<nmes2pts_known) out++;
    
    //check out
    if(out==nmes2pts_known)
      {
	master_fprintf(stderr,"Erorr, unkwnown tag %s, use one in this list:\n",tag);
	for(int i=0;i<nmes2pts_known;i++) master_fprintf(stderr," %s\n",mes2pts_tag[i]);
	crash("See previous message");
      }
    
    return (mes2pts_known_t)out;
  }
  
  //read the list of mesons in terms of quarks
  void read_mes2pts_contr_pars()
  {
    int nmes2pts_contr;
    read_str_int("NMes2PtsContr",&nmes2pts_contr);
    for(int i=0;i<nmes2pts_contr;i++)
      {
	char name[1024];
	read_str(name,1024);
	char q_name[2][1024];
	for(int iq=0;iq<2;iq++)
	  {
	    read_str(q_name[iq],1024);
	    if(Q.find(q_name[iq])==Q.end()) crash("unable to find q%d %s",iq,q_name[iq]);
	  }
	mes2pts_contr_map.push_back(mes_contr_map_t(name,q_name[0],q_name[1]));
      }
    
    if(nmes2pts_contr) read_mes2pts_contr_gamma_list();
  }
  
  //read the list of meson contraction asked
  void read_mes2pts_contr_gamma_list()
  {
    int nmes_gamma_contr;
    read_str_int("NGammaContr",&nmes_gamma_contr);
    for(int i=0;i<nmes_gamma_contr;i++)
      {
	switch(read_2pts_tag())
	  {
	  case P5P5: mes_gamma_list.push_back(idirac_pair_t(5,5));                                 break;
	  case P5GI: for(int ig=0;ig<16;ig++) mes_gamma_list.push_back(idirac_pair_t(5,ig));       break;
	  case GIP5: for(int ig=0;ig<16;ig++) mes_gamma_list.push_back(idirac_pair_t(ig,5));       break;
	  case S0GI: for(int ig=0;ig<16;ig++) mes_gamma_list.push_back(idirac_pair_t(0,ig));       break;
	  case GIS0: for(int ig=0;ig<16;ig++) mes_gamma_list.push_back(idirac_pair_t(ig,0));       break;
	  case V0V0: mes_gamma_list.push_back(idirac_pair_t(4,4));                                 break;
	  case AKAK: for(int mu=1;mu<=3;mu++) mes_gamma_list.push_back(idirac_pair_t(mu+5,mu+5));  break;
	  case VKVK: for(int mu=1;mu<=3;mu++) mes_gamma_list.push_back(idirac_pair_t(mu,mu));      break;
	  case VKTK: for(int mu=1;mu<=3;mu++) mes_gamma_list.push_back(idirac_pair_t(mu,mu+9));    break;
	  case TKVK: for(int mu=1;mu<=3;mu++) mes_gamma_list.push_back(idirac_pair_t(mu+9,mu));    break;
	  case TKTK: for(int ig=10;ig<=12;ig++) mes_gamma_list.push_back(idirac_pair_t(ig,ig));    break;
	  case BKBK: for(int ig=13;ig<=15;ig++) mes_gamma_list.push_back(idirac_pair_t(ig,ig));    break;
	  case V0P5: mes_gamma_list.push_back(idirac_pair_t(4,5));                                 break;
	  case VKP5: for(int mu=1;mu<=3;mu++) mes_gamma_list.push_back(idirac_pair_t(mu,5));       break;
	  case S0S0: mes_gamma_list.push_back(idirac_pair_t(0,0));                                 break;
	  case A0A0: mes_gamma_list.push_back(idirac_pair_t(9,9));                                 break;
	  case AKBK: for(int mu=1;mu<=3;mu++) mes_gamma_list.push_back(idirac_pair_t(mu+5,mu+12)); break;
	  case BKAK: for(int mu=1;mu<=3;mu++) mes_gamma_list.push_back(idirac_pair_t(mu+12,mu+5)); break;
	  case S0V0: mes_gamma_list.push_back(idirac_pair_t(0,4));                                 break;
	  case V0S0: mes_gamma_list.push_back(idirac_pair_t(4,0));                                 break;
	  default: crash("unknown meson_contr");
	  }
      }
  }
  
  //read the list of baryons in terms of quarks
  void read_bar2pts_contr_pars()
  {
    int nbar2pts_contr;
    read_str_int("NBar2PtsContr",&nbar2pts_contr);
    if(nbar2pts_contr and (!diluted_col_source or !diluted_spi_source)) crash("to make baryon contractions you need diluted color and spin");
    if(nbar2pts_contr)
      {
	read_str_int("ComputeOctet",&compute_octet);
	read_str_int("ComputeDecuplet",&compute_decuplet);
      }
    for(int i=0;i<nbar2pts_contr;i++)
      {
	char name[1024];
	read_str(name,1024);
	char q_name[3][1024];
	for(int iq=0;iq<3;iq++)
	  {
	    read_str(q_name[iq],1024);
	    if(Q.find(q_name[iq])==Q.end()) crash("unable to find q%d %s",iq,q_name[iq]);
	  }
	bar2pts_contr_map.push_back(bar_triplet_t(name,q_name[0],q_name[1],q_name[2]));
      }
  }
  
  //read the list of "handcuffs" contraction
  void read_handcuffs_contr_pars()
  {
    int nhand_contr;
    read_str_int("NHandcuffsContr",&nhand_contr);
    
    if(nhand_contr)
      {
	int nhand_sides;
	read_str_int("NHandcuffsSides",&nhand_sides);
	
	//read the sides
	for(int i=0;i<nhand_sides;i++)
	  {
	    char name[1024];
	    read_str(name,1024);
	    int igamma;
	    read_int(&igamma);
	    char bw[1024];
	    char fw[1024];
	    read_str(bw,1024);
	    read_str(fw,1024);
	    int store;
	    read_int(&store);
	    
	    handcuffs_side_map.push_back(handcuffs_side_map_t(name,igamma,bw,fw,store));
	  }
	
	//read the sides combo
	for(int i=0;i<nhand_contr;i++)
	  {
	    char name[1024];
	    read_str(name,1024);
	    char left[1024];
	    char right[1024];
	    read_str(left,1024);
	    read_str(right,1024);
	    
	    handcuffs_map.push_back(handcuffs_map_t(name,left,right));
	  }
      }
  }
  
  //read the parameters of the propagators to Fourier transform
  void read_fft_prop_pars()
  {
    //read the number of props to fft
    int nfft_props;
    read_str_int("NFftProps",&nfft_props);
    
    if(nfft_props)
      {
	//read the list of propagators
	for(int ifft_prop=0;ifft_prop<nfft_props;ifft_prop++)
	  {
	    char tag[1024];
	    read_str(tag,1024);
	    
	    if(Q.find(tag)==Q.end()) crash("unable to find %s",tag);
	    fft_prop_list.push_back(tag);
	  }
	
	//read the number of ranges
	int nfft_ranges;
	read_str_int("NFftRanges",&nfft_ranges);
	
	if(nfft_ranges>0)
	  {
	    std::vector<std::pair<fft_mom_range_t,double>> fft_mom_range_list;
	    
	    //read all the ranges
	    for(int irange=0;irange<nfft_ranges;irange++)
	      {
		fft_mom_range_t fft_mom_range;
		int L[2],T[2];
		read_str_int("L",&L[0]);
		read_int(&L[1]);
		read_str_int("T",&T[0]);
		read_int(&T[1]);
		
		//init the offset and width from range interval
		fft_mom_range.offs[0]=T[0];
		fft_mom_range.width[0]=T[1]-T[0]+1;
		for(int i=1;i<NDIM;i++)
		  {
		    fft_mom_range.offs[i]=L[0];
		    fft_mom_range.width[i]=L[1]-L[0]+1;
		  }
		
		//read the democratic filter
		double p4_fr_p22_max;
		read_str_double("P4FrP22Max",&p4_fr_p22_max);
		
		fft_mom_range_list.push_back(std::make_pair(fft_mom_range,p4_fr_p22_max));
	      }
	    
	    //initialize the filters
	    init_fft_filter_from_range(fft_mom_range_list);
	  }
	else
	  for(int i=0;i<std::abs(nfft_ranges);i++)
	    {
	      char path[1024];
	      read_str_str("MomentaListPath",path,1024);
	      char suffix[1024];
	      read_str_str("Suffix",suffix,1024);
	      
	      init_fft_filterer_from_file(suffix,path);
	    }
      }
  }
}
