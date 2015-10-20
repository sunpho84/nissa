#include <nissa.hpp>

using namespace nissa;

//convention on gospel
const int follow_chris=0,follow_nazario=1;

//define types of quark propagator used
const int nins_kind=5;
enum insertion_t{                    ORIGINAL,  SCALAR,  PSEUDO,  PHOTON,   TADPOLE};
const char ins_name[nins_kind][20]={"ORIGINAL","SCALAR","PSEUDO","PHOTON", "TADPOLE"};
const int nqprop_kind=6;
enum qprop_t{                           PROP_0,  PROP_S,  PROP_P,  PROP_T,  PROP_PHOTON,  PROP_PHOTON2};
const char prop_name[nqprop_kind][20]={"PROP_0","PROP_S","PROP_P","PROP_T","PROP_PHOTON","PROP_PHOTON2"};

//map the source, the destination and the insertion for each propagator
const qprop_t prop_map[nqprop_kind]=         {PROP_0,   PROP_S, PROP_P, PROP_T,   PROP_PHOTON,  PROP_PHOTON2};
const insertion_t insertion_map[nqprop_kind]={ORIGINAL, SCALAR, PSEUDO, TADPOLE,  PHOTON,       PHOTON};
const qprop_t source_map[nqprop_kind]=       {PROP_0,   PROP_0, PROP_0, PROP_0,   PROP_0,       PROP_PHOTON};
const char prop_abbr[]=                       "0"       "S"     "P"     "T"       "L"           "M";

//sign of the lepton momentum
const int norie=2;
const int sign_orie[2]={-1,+1};

const int nins=2;
const int nrev=2;
