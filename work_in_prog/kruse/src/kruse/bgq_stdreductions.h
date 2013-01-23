/*
 * bgq_stdreductions.h
 *
 *  Created on: Nov 21, 2012
 *      Author: meinersbur
 */

#ifndef BGQ_STDREDUCTIONS_H_
#define BGQ_STDREDUCTIONS_H_

#include "bgq_spinorfield.h"


complex_double bgq_spinorfield_innerprod_local(bool isOdd, bgq_weylfield_controlblock *field1, bgq_weylfield_controlblock *field2);
complex_double bgq_spinorfield_innerprod_global(bool isOdd, bgq_weylfield_controlblock *field1, bgq_weylfield_controlblock *field2);

double bgq_spinorfield_innerprod_r_local(bool isOdd, bgq_weylfield_controlblock *field1, bgq_weylfield_controlblock *field2);
double bgq_spinorfield_innerprod_r_global(bool isOdd, bgq_weylfield_controlblock *field1, bgq_weylfield_controlblock *field2);

double bgq_spinorfield_sqrnorm_local(tristate isOdd, bgq_weylfield_controlblock *field);
double bgq_spinorfield_sqrnorm_global(tristate isOdd, bgq_weylfield_controlblock *field);

#endif /* BGQ_STDREDUCTIONS_H_ */
