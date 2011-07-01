#define real double
#define complexx compled
#define sqrtx sqrt
#define rsqrtx rsqrt

su3g in;
sploat_to_double_vec((double*)in,in_high,in_low,8);

#include "gauge_decompression.c"

#undef real
#undef complexx
#undef sqrtx
#undef rsqrtx
