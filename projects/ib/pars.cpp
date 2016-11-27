#include <nissa.hpp>

#define EXTERN_PARS
#include "pars.hpp"

#include "prop.hpp"
#include "conf.hpp"
#include "contr.hpp"

namespace nissa
{
  //read all the parameters to contract with leptons
  void read_meslep_contr_pars()
  {
    //Leptons
    if(twisted_run) read_str_int("NMesLepQ1Q2LepmassMesMass",&nquark_lep_combos);
    else            read_str_int("NMesLepQ1Q2LepkappaMesMass",&nquark_lep_combos);
    
    if(nquark_lep_combos)  read_gospel_convention();
    
    lep_contr_iq1=nissa_malloc("lep_contr_iq1",nquark_lep_combos,int);
    lep_contr_iq2=nissa_malloc("lep_contr_iq2",nquark_lep_combos,int);
    leps=nissa_malloc("leps",nquark_lep_combos,tm_quark_info);
    lep_energy=nissa_malloc("lep_energy",nquark_lep_combos,double);
    neu_energy=nissa_malloc("neu_energy",nquark_lep_combos,double);
    for(int il=0;il<nquark_lep_combos;il++)
      {
	//read quarks identfiying the mesons
	read_int(lep_contr_iq1+il);
	read_int(lep_contr_iq2+il);
	
	//if not pure wilson read mass
	if(!twisted_run) leps[il].mass=0;
	else             read_double(&leps[il].mass);
	
	//antiperiodic or periodic
	leps[il].bc[0]=QUARK_BOUND_COND;
	
	//maximal twist (if tm), otherwise read kappa
	if(!twisted_run) read_double(&leps[il].kappa);
	else             leps[il].kappa=0.125;
	leps[il].r=0;
	
	//read the mass of the meson (that must have been determined outside)
	double mes_mass;
	read_double(&mes_mass);
	
	//set initial value of bc and check kinematic
	for(int i=1;i<NDIM;i++) leps[il].bc[i]=0;
	if(tm_quark_energy(leps[il],0)>=mes_mass) crash("initial state is lighter (%lg) than final state at rest (%lg)!",mes_mass,tm_quark_energy(leps[il],0));
	
	//compute meson momentum and bc
	double err;
	master_printf("Resolving kinematical condition for combination of quarks %d/%d\n",il+1,nquark_lep_combos);
	do
	  {
	    //compute the error
	    double lep_energy=tm_quark_energy(leps[il],0);
	    double neu_energy=naive_massless_quark_energy(leps[il].bc,0);
	    err=lep_energy+neu_energy-mes_mass;
	    //compute the derivative
	    double eps=1e-8;
	    for(int i=1;i<NDIM;i++) leps[il].bc[i]+=eps;
	    double der=(tm_quark_energy(leps[il],0)+naive_massless_quark_energy(leps[il].bc,0)-mes_mass-err)/eps;
	    for(int i=1;i<NDIM;i++) leps[il].bc[i]-=eps+err/der;
	    
	    master_printf("  lep_e: %+010.10lg, neu_e: %+010.10lg, mes_mass: %lg, error: %lg, der: %lg\n",lep_energy,neu_energy,mes_mass,err,der);
	  }
	while(fabs(err)>1e-14);
	
	//write down energy
	lep_energy[il]=tm_quark_energy(leps[il],0);
	neu_energy[il]=naive_massless_quark_energy(leps[il].bc,0);
	master_printf(" ilepton %d, lepton energy: %lg, neutrino energy: %lg\n",il,lep_energy[il],neu_energy[il]);
	master_printf(" lep+neut energy: %lg\n",lep_energy[il]+neu_energy[il]);
	master_printf(" bc: %+016.016lg\n\n",leps[il].bc[1]);
      }
  }
  
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
    //master_printf("Tadpole: {%lg, %lg, %lg, %lg}\n",tadpole[0],tadpole[1],tadpole[2],tadpole[3]);
  }
  
  //meson tags
  const int nmes2pts_known=9;
  enum mes2pts_known_t                      { P5P5 , P5GI , V0V0 , AKAK , VKVK , VKTK , TKVK , TKTK , BKBK };
  const char mes2pts_tag[nmes2pts_known][5]={"P5P5","P5GI","V0V0","AKAK","VKVK","VKTK","TKVK","TKTK","BKBK"};
  mes2pts_known_t read_2pts_tag()
  {
    //read the tag
    char tag[10];
    read_str(tag,10);
    
    //convert to int
    int out=0;
    while(strcasecmp(tag,mes2pts_tag[out]) && out<nmes2pts_known) out++;
    
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
	  case P5P5: mes_gamma_list.push_back(idirac_pair_t(5,5));                                break;
	  case P5GI: for(int ig=0;ig<16;ig++) mes_gamma_list.push_back(idirac_pair_t(5,ig));      break;
	  case V0V0: mes_gamma_list.push_back(idirac_pair_t(4,4));                                break;
	  case AKAK: for(int mu=1;mu<=3;mu++) mes_gamma_list.push_back(idirac_pair_t(mu+5,mu+5)); break;
	  case VKVK: for(int mu=1;mu<=3;mu++) mes_gamma_list.push_back(idirac_pair_t(mu,mu));     break;
	  case VKTK: for(int mu=1;mu<=3;mu++) mes_gamma_list.push_back(idirac_pair_t(mu,mu+9));   break;
	  case TKVK: for(int mu=1;mu<=3;mu++) mes_gamma_list.push_back(idirac_pair_t(mu+9,mu));   break;
	  case TKTK: for(int ig=10;ig<=12;ig++) mes_gamma_list.push_back(idirac_pair_t(ig,ig));   break;
	  case BKBK: for(int ig=13;ig<=15;ig++) mes_gamma_list.push_back(idirac_pair_t(ig,ig));   break;
	  default: crash("unknown meson_contr");
	  }
      }
  }
  
  //read the list of baryons in terms of quarks
  void read_bar2pts_contr_pars()
  {
    int nbar2pts_contr;
    read_str_int("NBar2PtsContr",&nbar2pts_contr);
    if(nbar2pts_contr&&((!diluted_col_source)||(!diluted_spi_source))) crash("to make baryon contractions you need diluted color and spin");
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
}
