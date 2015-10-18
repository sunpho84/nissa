#include <nissa.hpp>

using namespace nissa;

//convention on gospel
const int follow_chris=0,follow_nazario=1;

//kind of photon sources
const int nphi_eta_alt=3;
const int iphi=0,ieta=1,ialt=2;
const char photon_field_name[nphi_eta_alt][4]={"phi","eta","alt"};

//define types of quark propagator used
const int nins_kind=7;
enum insertion_t{                    ORIGINAL,  SCALAR,  PSEUDO,  STOCH_PHI,  STOCH_ETA,  STOCH_ALT,   TADPOLE};
const char ins_name[nins_kind][20]={"ORIGINAL","SCALAR","PSEUDO","STOCH_PHI","STOCH_ETA","STOCH_ALT", "TADPOLE"};
const int nqprop_kind=9;
enum qprop_t{                           PROP_0,  PROP_S,  PROP_P,  PROP_PHI,  PROP_ETA,  PROP_PHIETA,  PROP_T,  PROP_ALT,  PROP_ALT_ALT};
const char prop_name[nqprop_kind][20]={"PROP_0","PROP_S","PROP_P","PROP_PHI","PROP_ETA","PROP_PHIETA","PROP_T","PROP_ALT","PROP_ALT_ALT"};
const qprop_t PROP_PHI_ETA_ALT[nphi_eta_alt]={PROP_PHI,PROP_ETA,PROP_ALT};

//map the source, the destination and the insertion for each propagator
const qprop_t prop_map[nqprop_kind]=         {PROP_0,   PROP_S, PROP_P, PROP_PHI,  PROP_ETA,  PROP_PHIETA, PROP_T,   PROP_ALT,  PROP_ALT_ALT};
const insertion_t insertion_map[nqprop_kind]={ORIGINAL, SCALAR, PSEUDO, STOCH_PHI, STOCH_ETA, STOCH_ETA,   TADPOLE,  STOCH_ALT, STOCH_ALT};
const qprop_t source_map[nqprop_kind]=       {PROP_0,   PROP_0, PROP_0, PROP_0,    PROP_0,    PROP_PHI,    PROP_0,   PROP_0,    PROP_ALT};
const char prop_abbr[]=                       "0"       "S"     "P"     "A"        "B"        "X"          "T"       "L"        "M";

//compute the eigenvalues of (1-+g0)/2
const double W=1/sqrt(2);
const spin ompg0_eig[2][2]={{{{+W, 0},{ 0, 0},{+W, 0},{ 0, 0}},
		       {{ 0, 0},{+W, 0},{ 0, 0},{+W, 0}}},
		      {{{+W, 0},{ 0, 0},{-W, 0},{ 0, 0}},
		       {{ 0, 0},{+W, 0},{ 0, 0},{-W, 0}}}};

//sign of the lepton momentum
const int norie=2;
const int sign_orie[2]={-1,+1};

const int nins=2;
const int nrev=2;
