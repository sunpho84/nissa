#ifndef BGQ_COMM_H_
#define BGQ_COMM_H_

#include "bgq_spinorfield.h"
#include "bgq_field.h"
#include "bgq_utils.h"

#include <stdint.h>

#ifdef __cplusplus
extern "C"{
#endif

#ifndef BGQ_COMM_C_
#define EXTERN_INLINE EXTERN_INLINE_DECLARATION
#define EXTERN_FIELD extern
#define EXTERN_INIT(val)
#else
#define EXTERN_INLINE EXTERN_INLINE_DEFINITION
#define EXTERN_FIELD
#define EXTERN_INIT(val) = (val)
#endif


void bgq_comm_mpi_init(void);
void bgq_comm_spi_init(void);

//TODO: inline?
void bgq_comm_recv(bool nospi, bool sloppy, bgq_weylfield_controlblock *targetfield);
void bgq_comm_send(void);
void bgq_comm_wait(void);





EXTERN_FIELD uint8_t *g_bgq_sec_comm;
EXTERN_FIELD uint8_t *g_bgq_sec_comm_float;
EXTERN_FIELD bgq_weyl_vec_double *g_bgq_sec_recv_double[PHYSICAL_LD];
EXTERN_FIELD bgq_weyl_vec_double *g_bgq_sec_send_double[PHYSICAL_LD];

EXTERN_FIELD bgq_weyl_vec_double *g_bgq_sec_temp_tup_double;
EXTERN_FIELD bgq_weyl_vec_double *g_bgq_sec_temp_tdown_double;


EXTERN_FIELD bgq_weyl_vec_float *g_bgq_sec_recv_float[PHYSICAL_LD];
EXTERN_FIELD bgq_weyl_vec_float *g_bgq_sec_send_float[PHYSICAL_LD];

EXTERN_FIELD bgq_weyl_vec_float *g_bgq_sec_temp_tup_float;
EXTERN_FIELD bgq_weyl_vec_float *g_bgq_sec_temp_tdown_float;


#define g_bgq_sec_recv NAME2(g_bgq_sec_recv,PRECISION)
#define g_bgq_sec_send NAME2(g_bgq_sec_send,PRECISION)

#define g_bgq_sec_temp_tup NAME2(g_bgq_sec_temp_tup,PRECISION)
#define g_bgq_sec_temp_tdown NAME2(g_bgq_sec_temp_tdown,PRECISION)


#undef EXTERN_INLINE
#undef EXTERN_FIELD
#undef EXTERN_INIT

#ifdef __cplusplus
}
#endif


#endif /* BGQ_COMM_H_ */
