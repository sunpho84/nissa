/*
 * bgq_utils.h
 *
 *  Created on: Aug 4, 2012
 *      Author: meinersbur
 */


#define BGQ_UTILS_C_
#include "bgq_utils.h"

#include <mpi.h>

#include <stdlib.h>
#include <stdio.h>
#include <string.h>


void *malloc_aligned(size_t size, size_t alignment) {
	void *result = NULL;
	if(size==0) size=1;
	int errcode = posix_memalign(&result, alignment, size);
	if (errcode != 0) {
		fprintf(stderr, "malloc returned %d\n", errcode);
		exit(10);
	}
	//memset(result, 0, size);
	return result;
}

char *x = NULL;
void opaque_func_call() {
	// Do some no-op the compiler cannot opt away
	if (x && *x) {
		*x = '\0';
	}

	return;
	// clear the caches
	char *p = malloc(20*1024*1024);
	char *end = p+20*1024*1024;
	for (char *q = p; q<end; q+=1) {
		*q='\0';
	}
	free(p);
}
