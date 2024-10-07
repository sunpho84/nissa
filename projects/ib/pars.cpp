#include <nissa.hpp>

#define EXTERN_PARS
#include "pars.hpp"

#include "prop.hpp"
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
    if(strncasecmp(photon_gauge_str,"FEYNMAN",100)==0) photon.which_gauge=gauge_info::FEYNMAN;
    else
      if(strncasecmp(photon_gauge_str,"LANDAU",100)==0) photon.which_gauge=gauge_info::LANDAU;
      else
	if(strncasecmp(photon_gauge_str,"COULOMB",100)==0) photon.which_gauge=gauge_info::COULOMB;
	else crash("Unkwnown photon gauge: %s",photon_gauge_str);
    
    //discretization for photon propagator
    char photon_discrete_str[100];
    read_str_str("PhotonDiscretization",photon_discrete_str,100);
    if(strncasecmp(photon_discrete_str,"WILSON",100)==0) photon.c1=WILSON_C1;
    else
      if(strncasecmp(photon_discrete_str,"TLSYM",100)==0) photon.c1=TLSYM_C1;
      else crash("Unkwnown photon discretization: %s",photon_discrete_str);
    
    //compute the tadpole summing all momentum
    tadpole=compute_tadpole(photon);
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
	  read_str(q_name[iq],1024);
	
	for(int icopy=0;icopy<ncopies;icopy++)
	  {
	    char suffix[128]="";
	    if(ncopies>1) sprintf(suffix,"_copy%d",icopy);
	    
	    char full_name[1024+129];
	    sprintf(full_name,"%s%s",name,suffix);
	    
	    char q_full_name[2][1024+129];
	    for(int iq=0;iq<2;iq++)
	      sprintf(q_full_name[iq],"%s%s",q_name[iq],suffix);
	    for(int iq=0;iq<2;iq++)
	      if(Q.find(q_full_name[iq])==Q.end()) crash("unable to find q%d %s",iq,q_full_name[iq]);
	    
	    mes2pts_contr_map.push_back(mes_contr_map_t(full_name,q_full_name[0],q_full_name[1]));
	    for(int iq=0;iq<2;iq++)
	      propsNeededToContr.insert(q_full_name[iq]);
	  }
      }
    
    if(nmes2pts_contr) read_mes2pts_contr_gamma_list();
  }
  
  /// List of gamma used for source or sink
  struct LocBilinear
  {
    /// Covnerted list
    const std::vector<int> list;
    
    /// Letter to identify
    const char letter;
    
    /// Convert to int a char
    static int CliffOfChar(const char& g)
    {
      static constexpr char CliffMap[8]=
	"SVPATBG";
      
      for(int i=0;i<7;i++)
	if(CliffMap[i]==g)
	  return
	    i;
      
      crash("Cannot convert gamma: %c",g);
      
      return
	{};
    }
    
  static std::vector<int> getList(const char& g,
				  const char& letter)
    {
      static const std::vector<int> literal[7]={
	{0},
	{1,2,3},
	{5},
	{6,7,8},
	{10,11,12},
	{13,14,15},
	{0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15}};
      
      constexpr static bool ignoreLetter[7]=
		  {1,0,1,0,0,0,1};
      
      const int c=
	CliffOfChar(g);
      
      if(ignoreLetter[c] or not isdigit(letter))
	return
	  literal[c];
      
      static const int minNumeric[7]={0,0,0,0,1,1,0};
      
      const int m=minNumeric[c];
      
      const int i=
	letter-'0'-m;
      
      static const std::vector<int> numeric[7]={
	{0},
	{4,1,2,3},
	{5},
	{9,6,7,8},
	{10,11,12},
	{13,14,15},
	{0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15}};
      
      if(i<0 or i>=(int)numeric[c].size())
	crash("Error, letter %c converts to int %d not in range [%d:%d]",letter,i,m,numeric[c].size()+m);
      
      return
	{numeric[c][i]};
    }
  
  LocBilinear(const char& g,
	    const char& letter)
    : list(getList(g,letter)),letter(letter)
    {
    }
  };
  
  //read the list of meson contraction asked
  void read_mes2pts_contr_gamma_list()
  {
    int nmes_gamma_contr;
    read_str_int("NGammaContr",&nmes_gamma_contr);
    for(int i=0;i<nmes_gamma_contr;i++)
      {
	char tag[128];
	read_str(tag,128);
	
	if(strlen(tag)!=4)
	  crash("Error, string length %lu different from 4",strlen(tag));
	
	LocBilinear sink(tag[0],tag[1]);
	LocBilinear source(tag[2],tag[3]);
	
	if(source.letter!=sink.letter)
	  for(const int& iSink : sink.list)
	    for(const int& iSource : source.list)
	      mes_gamma_list.push_back({iSink,iSource});
	else
	  {
	    if(source.list.size()!=sink.list.size())
	      crash("Error, sizes of source %lu does not agree with size of sink %lu",source.list.size(),sink.list.size());
	  
	  for(size_t i=0;i<source.list.size();i++)
	    mes_gamma_list.push_back({sink.list[i],source.list[i]});
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
	  read_str(q_name[iq],1024);
	
	for(int icopy=0;icopy<ncopies;icopy++)
	  {
	    char suffix[128]="";
	    if(ncopies>1) sprintf(suffix,"_copy%d",icopy);
	    
	    char full_name[1024+129];
	    sprintf(full_name,"%s%s",name,suffix);
	    char q_full_name[3][1024+129];
	    for(int iq=0;iq<3;iq++)
	      sprintf(q_full_name[iq],"%s%s",q_name[iq],suffix);
	    for(int iq=0;iq<3;iq++)
	      if(Q.find(q_full_name[iq])==Q.end()) crash("unable to find q%d %s",iq,q_full_name[iq]);
	    bar2pts_contr_map.push_back(bar_triplet_t(full_name,q_full_name[0],q_full_name[1],q_full_name[2]));
	    for(int iq=0;iq<3;iq++)
	      propsNeededToContr.insert(q_full_name[iq]);
	  }
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
	    
	    for(int icopy=0;icopy<ncopies;icopy++)
	      {
		char suffix[128]="";
		if(ncopies>1) sprintf(suffix,"_copy%d",icopy);
		
		char full_name[1024+129];
		char bw_full[1024+129];
		char fw_full[1024+129];
		sprintf(full_name,"%s%s",name,suffix);
		sprintf(bw_full,"%s%s",bw,suffix);
		sprintf(fw_full,"%s%s",fw,suffix);
		if(Q.find(bw_full)==Q.end()) crash("for bubble \'%s\' the first propagator \'%s\' is not present",name,bw_full);
		if(Q.find(fw_full)==Q.end()) crash("for bubble \'%s\' the second propagator \'%s\' is not present",name,fw_full);
		handcuffs_side_map.push_back(handcuffs_side_map_t(full_name,igamma,bw_full,fw_full,store));
		for(auto& q : {bw_full,fw_full})
		  propsNeededToContr.insert(q);
	      }
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
	    
	    for(int icopy=0;icopy<ncopies;icopy++)
	      {
		char suffix[128]="";
		if(ncopies>1) sprintf(suffix,"_copy%d",icopy);
		
		char full_name[1024+129];
		char left_full[1024+129];
		char right_full[1024+129];
		sprintf(full_name,"%s%s",name,suffix);
		sprintf(left_full,"%s%s",left,suffix);
		sprintf(right_full,"%s%s",right,suffix);
		
		//check if left and right is present in handcuffs_size_map
		bool left_hand_found=false;
		bool right_hand_found=false;
		for(auto &hand : handcuffs_side_map)
		  {
		    if(hand.name==left_full) left_hand_found=true;
		    if(hand.name==right_full) right_hand_found=true;
		    if(left_hand_found && right_hand_found) break;
		  }
		
		if(!left_hand_found)   crash("for handcuffs \'%s\' the left bubble \'%s\' is not present",name,left_full);
		if(!right_hand_found)  crash("for handcuffs \'%s\' the right bubble \'%s\' is not present",name,right_full);
		
		handcuffs_map.push_back(handcuffs_map_t(full_name,left_full,right_full));
	      }
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
	    
	    for(int icopy=0;icopy<ncopies;icopy++)
	      {
		char suffix[128]="";
		if(ncopies>1) sprintf(suffix,"_copy%d",icopy);
		
		char tag_full[1024+129];
		sprintf(tag_full,"%s%s",tag,suffix);
		
		if(Q.find(tag_full)==Q.end()) crash("unable to find %s",tag_full);
		fft_prop_list.push_back(tag_full);
		propsNeededToContr.insert(tag_full);
	      }
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
