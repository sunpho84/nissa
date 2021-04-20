#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#define EXTERN_DIRAC
#include "dirac.hpp"

#include <math.h>
#include <string.h>

#include "complex.hpp"
#include "dirac.hpp"
#include "spin.hpp"
#include "su3.hpp"

#include "base/debug.hpp"
#include "geometry/geometry_lx.hpp"

namespace nissa
{
  //Initialize the gamma matrix base and the rotators
  // base_gamma[ 0] = identity
  // base_gamma[ 1] = gamma_1
  // base_gamma[ 2] = gamma_2
  // base_gamma[ 3] = gamma_3
  // base_gamma[ 4] = gamma_0
  // base_gamma[ 5] = gamma_5
  // base_gamma[ 6] = gamma_15
  // base_gamma[ 7] = gamma_25
  // base_gamma[ 8] = gamma_35
  // base_gamma[ 9] = gamma_05
  // base_gamma[10] = gamma_01
  // base_gamma[11] = gamma_02
  // base_gamma[12] = gamma_03
  // base_gamma[13] = gamma_23
  // base_gamma[14] = gamma_31
  // base_gamma[15] = gamma_12
  
  void init_base_gamma()
  {
    const double rad2=1./sqrt(2);
    
    init_dirac(base_gamma+ 0,  0,1,0  , 1,1,0  , 2,1,0  , 3,1,0 );
    init_dirac(base_gamma+ 1,  3,0,-1 , 2,0,-1 , 1,0,1  , 0,0,1 );
    init_dirac(base_gamma+ 2,  3,-1,0 , 2,1,0  , 1,1,0  , 0,-1,0);
    init_dirac(base_gamma+ 3,  2,0,-1 , 3,0,1  , 0,0,1  , 1,0,-1);
    init_dirac(base_gamma+ 4,  2,-1,0 , 3,-1,0 , 0,-1,0 , 1,-1,0);
    init_dirac(base_gamma+ 5,  0,1,0  , 1,1,0  , 2,-1,0 , 3,-1,0);
    init_dirac(base_gamma+ 6,  3,0,1  , 2,0,1  , 1,0,1  , 0,0,1 );
    init_dirac(base_gamma+ 7,  3,1,0  , 2,-1,0 , 1,1,0  , 0,-1,0);
    init_dirac(base_gamma+ 8,  2,0,1  , 3,0,-1 , 0,0,1  , 1,0,-1);
    init_dirac(base_gamma+ 9,  2,1,0  , 3,1,0  , 0,-1,0 , 1,-1,0);
    init_dirac(base_gamma+10,  1,0,-1 , 0,0,-1 , 3,0,1  , 2,0,1 );
    init_dirac(base_gamma+11,  1,-1,0 , 0,1,0  , 3,1,0  , 2,-1,0);
    init_dirac(base_gamma+12,  0,0,-1 , 1,0,1  , 2,0,1  , 3,0,-1);
    init_dirac(base_gamma+13,  1,0,1  , 0,0,1  , 3,0,1  , 2,0,1 );
    init_dirac(base_gamma+14,  1,1,0  , 0,-1,0 , 3,1,0  , 2,-1,0);
    init_dirac(base_gamma+15,  0,0,1  , 1,0,-1 , 2,0,1  , 3,0,-1);
    
    init_dirac(&Pplus ,          0,rad2,rad2  , 1,rad2,rad2  , 2,rad2,-rad2 , 3,rad2,-rad2);
    init_dirac(&Pminus,          0,rad2,-rad2 , 1,rad2,-rad2 , 2,rad2,rad2  , 3,rad2,rad2 );
    
    dirac_prod(base_gamma+16,base_gamma+4,base_gamma+6);
    dirac_prod(base_gamma+17,base_gamma+4,base_gamma+7);
    dirac_prod(base_gamma+18,base_gamma+4,base_gamma+8);
    
    //the sigma mu nu in an anti-simmetric tensor
    int lmunu[6]={10,11,12,15,14,13};
    int cmunu[6]={1,1,1,1,-1,1};
    for(int imunu=0;imunu<6;imunu++)
      for(int d1=0;d1<NDIRAC;d1++)
	{
	  smunu_entr[d1][imunu][0]=cmunu[imunu]*base_gamma[lmunu[imunu]].entr[d1][0];
	  smunu_entr[d1][imunu][1]=cmunu[imunu]*base_gamma[lmunu[imunu]].entr[d1][1];
	  
	  smunu_pos[d1][imunu]=base_gamma[lmunu[imunu]].pos[d1];
	}
    
    //1+gamma_mu, 1-gamma_mu
    for(int mu=0;mu<NDIRAC;mu++) //nasty
      {
	spinspin_put_to_id(opg[mu]);
	spinspin_put_to_id(omg[mu]);
	
	spinspin_dirac_summ_the_prod_double(opg[mu],&(base_gamma[igamma_of_mu(Dir(mu))]),+1);
	spinspin_dirac_summ_the_prod_double(omg[mu],&(base_gamma[igamma_of_mu(Dir(mu))]),-1);
      }
  }
}
