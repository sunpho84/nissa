/*
 * bgq_readWeyllayout.inc.c
 *
 *  Created on: Nov 13, 2012
 *      Author: meinersbur
 */



#ifndef BGQ_READWEYLLAYOUT_INC_
#include "bgq_qpx.h"
#include "bgq_spinorfield.h"

#define PRECISION double

void bgq_readWeyllayout(bgq_su3_spinor_params(/*out*/spinor), bgq_weylsite *weylsite, ucoord t1, ucoord t2, ucoord x, ucoord y, ucoord z)
#endif

#ifndef BGQ_READWEYLLAYOUT_INSERTPREFETCH
#define BGQ_READWEYLLAYOUT_INSERTPREFETCH
#endif

{

	bgq_su3_weylnext_prefetch(weylsite);

	// T+ /////////////////////////////////////////////////////////////////////////
	{
		bgq_su3_weyl_decl(weylnext_tup);
		bgq_qvlfuxa(weylnext_tup_v0_c0, weylsite, BGQ_VECTORSIZE);
		bgq_qvlfuxa(weylnext_tup_v0_c1, weylsite, BGQ_VECTORSIZE);
		bgq_qvlfuxa(weylnext_tup_v0_c2, weylsite, BGQ_VECTORSIZE);
		bgq_qvlfuxa(weylnext_tup_v1_c0, weylsite, BGQ_VECTORSIZE);
		bgq_qvlfuxa(weylnext_tup_v1_c1, weylsite, BGQ_VECTORSIZE);
		bgq_qvlfuxa(weylnext_tup_v1_c2, weylsite, BGQ_VECTORSIZE);
				bgq_weylqpx_expect(weylnext_tup, t1, t2, x, y, z, TUP, false);
				bgq_setdesc(BGQREF_TUP_RECV, "BGQREF_TUP_RECV");
				bgq_setbgqvalue(t1, x, y, z, BGQREF_TUP_RECV, bgq_cmplxval1(weylnext_tup_v0_c0));
				bgq_setbgqvalue(t2, x, y, z, BGQREF_TUP_RECV, bgq_cmplxval2(weylnext_tup_v0_c0));
		bgq_su3_expand_weyl_tup(spinor, weylnext_tup);
	}

	bgq_su3_weylnext_prefetch(weylsite);

	// T- /////////////////////////////////////////////////////////////////////////
	{
		bgq_su3_weyl_decl(weylnext_tdown);
		bgq_qvlfuxa(weylnext_tdown_v0_c0, weylsite, BGQ_VECTORSIZE);
		bgq_qvlfuxa(weylnext_tdown_v0_c1, weylsite, BGQ_VECTORSIZE);
		bgq_qvlfuxa(weylnext_tdown_v0_c2, weylsite, BGQ_VECTORSIZE);
		bgq_qvlfuxa(weylnext_tdown_v1_c0, weylsite, BGQ_VECTORSIZE);
		bgq_qvlfuxa(weylnext_tdown_v1_c1, weylsite, BGQ_VECTORSIZE);
		bgq_qvlfuxa(weylnext_tdown_v1_c2, weylsite, BGQ_VECTORSIZE);
				bgq_weylqpx_expect(weylnext_tdown, t1, t2, x, y, z, TDOWN, false);
				bgq_setdesc(BGQREF_TDOWN_RECV, "BGQREF_TDOWN_RECV");
				bgq_setbgqvalue(t1, x, y, z, BGQREF_TDOWN_RECV, bgq_cmplxval1(weylnext_tdown_v0_c0));
				bgq_setbgqvalue(t2, x, y, z, BGQREF_TDOWN_RECV, bgq_cmplxval2(weylnext_tdown_v0_c0));
		bgq_su3_accum_weyl_tdown(spinor, weylnext_tdown);
				bgq_setdesc(BGQREF_TDOWN_ACCUM, "BGQREF_TDOWN_ACCUM");
				bgq_setbgqvalue(t1, x, y, z, BGQREF_TDOWN_ACCUM, bgq_cmplxval1(spinor_v0_c0));
				bgq_setbgqvalue(t2, x, y, z, BGQREF_TDOWN_ACCUM, bgq_cmplxval2(spinor_v0_c0));
	}

	bgq_su3_weylnext_prefetch(weylsite);

	// X+ /////////////////////////////////////////////////////////////////////////
	{
		bgq_su3_weyl_decl(weylnext_xup);
		bgq_qvlfuxa(weylnext_xup_v0_c0, weylsite, BGQ_VECTORSIZE);
		bgq_qvlfuxa(weylnext_xup_v0_c1, weylsite, BGQ_VECTORSIZE);
		bgq_qvlfuxa(weylnext_xup_v0_c2, weylsite, BGQ_VECTORSIZE);
		bgq_qvlfuxa(weylnext_xup_v1_c0, weylsite, BGQ_VECTORSIZE);
		bgq_qvlfuxa(weylnext_xup_v1_c1, weylsite, BGQ_VECTORSIZE);
		bgq_qvlfuxa(weylnext_xup_v1_c2, weylsite, BGQ_VECTORSIZE);
				bgq_weylqpx_expect(weylnext_xup, t1, t2, x, y, z, XUP, false);
				bgq_setdesc(BGQREF_XUP_RECV, "BGQREF_XUP_RECV");
				bgq_setbgqvalue(t1, x, y, z, BGQREF_XUP_RECV, bgq_cmplxval1(weylnext_xup_v0_c0));
				bgq_setbgqvalue(t2, x, y, z, BGQREF_XUP_RECV, bgq_cmplxval2(weylnext_xup_v0_c0));
		bgq_su3_accum_weyl_xup(spinor, weylnext_xup);
				bgq_setdesc(BGQREF_XUP_ACCUM, "BGQREF_XDOWN_ACCUM");
				bgq_setbgqvalue(t1, x, y, z, BGQREF_XUP_ACCUM, bgq_cmplxval1(spinor_v0_c0));
				bgq_setbgqvalue(t2, x, y, z, BGQREF_XUP_ACCUM, bgq_cmplxval2(spinor_v0_c0));
	}

	bgq_su3_weylnext_prefetch(weylsite);

	// X- /////////////////////////////////////////////////////////////////////////
	{
		bgq_su3_weyl_decl(weylnext_xdown);
		bgq_qvlfuxa(weylnext_xdown_v0_c0, weylsite, BGQ_VECTORSIZE);
		bgq_qvlfuxa(weylnext_xdown_v0_c1, weylsite, BGQ_VECTORSIZE);
		bgq_qvlfuxa(weylnext_xdown_v0_c2, weylsite, BGQ_VECTORSIZE);
		bgq_qvlfuxa(weylnext_xdown_v1_c0, weylsite, BGQ_VECTORSIZE);
		bgq_qvlfuxa(weylnext_xdown_v1_c1, weylsite, BGQ_VECTORSIZE);
		bgq_qvlfuxa(weylnext_xdown_v1_c2, weylsite, BGQ_VECTORSIZE);
				bgq_weylqpx_expect(weylnext_xdown, t1, t2, x, y, z, XDOWN, false);
				bgq_setdesc(BGQREF_XDOWN_RECV, "BGQREF_XDOWN_RECV");
				bgq_setbgqvalue(t1, x, y, z, BGQREF_XDOWN_RECV, bgq_cmplxval1(weylnext_xdown_v0_c0));
				bgq_setbgqvalue(t2, x, y, z, BGQREF_XDOWN_RECV, bgq_cmplxval2(weylnext_xdown_v0_c0));
		bgq_su3_accum_weyl_xdown(spinor, weylnext_xdown);
				bgq_setdesc(BGQREF_XDOWN_ACCUM, "BGQREF_XDOWN_ACCUM");
				bgq_setbgqvalue(t1, x, y, z, BGQREF_XDOWN_ACCUM, bgq_cmplxval1(spinor_v0_c0));
				bgq_setbgqvalue(t2, x, y, z, BGQREF_XDOWN_ACCUM, bgq_cmplxval2(spinor_v0_c0));
	}

	bgq_su3_weylnext_prefetch(weylsite);

	// Y+ /////////////////////////////////////////////////////////////////////////
	{
		bgq_su3_weyl_decl(weylnext_yup);
		bgq_qvlfuxa(weylnext_yup_v0_c0, weylsite, BGQ_VECTORSIZE);
		bgq_qvlfuxa(weylnext_yup_v0_c1, weylsite, BGQ_VECTORSIZE);
		bgq_qvlfuxa(weylnext_yup_v0_c2, weylsite, BGQ_VECTORSIZE);
		bgq_qvlfuxa(weylnext_yup_v1_c0, weylsite, BGQ_VECTORSIZE);
		bgq_qvlfuxa(weylnext_yup_v1_c1, weylsite, BGQ_VECTORSIZE);
		bgq_qvlfuxa(weylnext_yup_v1_c2, weylsite, BGQ_VECTORSIZE);
				bgq_weylqpx_expect(weylnext_yup, t1, t2, x, y, z, YUP, false);
		bgq_su3_accum_weyl_yup(spinor, weylnext_yup);
	}

	bgq_su3_weylnext_prefetch(weylsite);

	// Y- /////////////////////////////////////////////////////////////////////////
	{
		bgq_su3_weyl_decl(weylnext_ydown);
		bgq_qvlfuxa(weylnext_ydown_v0_c0, weylsite, BGQ_VECTORSIZE);
		bgq_qvlfuxa(weylnext_ydown_v0_c1, weylsite, BGQ_VECTORSIZE);
		bgq_qvlfuxa(weylnext_ydown_v0_c2, weylsite, BGQ_VECTORSIZE);
		bgq_qvlfuxa(weylnext_ydown_v1_c0, weylsite, BGQ_VECTORSIZE);
		bgq_qvlfuxa(weylnext_ydown_v1_c1, weylsite, BGQ_VECTORSIZE);
		bgq_qvlfuxa(weylnext_ydown_v1_c2, weylsite, BGQ_VECTORSIZE);
				bgq_weylqpx_expect(weylnext_ydown, t1, t2, x, y, z, YDOWN, false);
		bgq_su3_accum_weyl_ydown(spinor, weylnext_ydown);
	}

	bgq_su3_weylnext_prefetch(weylsite);

	// Z+ /////////////////////////////////////////////////////////////////////////
	{
		bgq_su3_weyl_decl(weylnext_zup);
		bgq_qvlfuxa(weylnext_zup_v0_c0, weylsite, BGQ_VECTORSIZE);
		bgq_qvlfuxa(weylnext_zup_v0_c1, weylsite, BGQ_VECTORSIZE);
		bgq_qvlfuxa(weylnext_zup_v0_c2, weylsite, BGQ_VECTORSIZE);
		bgq_qvlfuxa(weylnext_zup_v1_c0, weylsite, BGQ_VECTORSIZE);
		bgq_qvlfuxa(weylnext_zup_v1_c1, weylsite, BGQ_VECTORSIZE);
		bgq_qvlfuxa(weylnext_zup_v1_c2, weylsite, BGQ_VECTORSIZE);
				bgq_weylqpx_expect(weylnext_zup, t1, t2, x, y, z, ZUP, false);
		bgq_su3_accum_weyl_zup(spinor, weylnext_zup);
	}

	BGQ_READWEYLLAYOUT_INSERTPREFETCH

	// Z- /////////////////////////////////////////////////////////////////////////
	{
		bgq_su3_weyl_decl(weylnext_zdown);
		bgq_qvlfuxa(weylnext_zdown_v0_c0, weylsite, BGQ_VECTORSIZE);
		bgq_qvlfuxa(weylnext_zdown_v0_c1, weylsite, BGQ_VECTORSIZE);
		bgq_qvlfuxa(weylnext_zdown_v0_c2, weylsite, BGQ_VECTORSIZE);
		bgq_qvlfuxa(weylnext_zdown_v1_c0, weylsite, BGQ_VECTORSIZE);
		bgq_qvlfuxa(weylnext_zdown_v1_c1, weylsite, BGQ_VECTORSIZE);
		bgq_qvlfuxa(weylnext_zdown_v1_c2, weylsite, BGQ_VECTORSIZE);
				bgq_weylqpx_expect(weylnext_zdown, t1, t2, x, y, z, ZDOWN, false);
		bgq_su3_accum_weyl_zdown(spinor, weylnext_zdown);
	}
}


#undef BGQ_READWEYLLAYOUT_INSERTPREFETCH
#undef BGQ_READWEYLLAYOUT_INC_
