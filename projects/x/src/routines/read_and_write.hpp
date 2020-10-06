#ifndef _READ_AND_WRITE_H
#define _READ_AND_WRITE_H
#include "../types/types.hpp"

void write_corr16(const char *path,corr16 *v,int prec);
void read_corr16(corr16 *v,const char *path);

#endif
