#ifndef _OPTIMIZED_MESONS_2PTS_H
#define _OPTIMIZED_MESONS_2PTS_H

void print_optimized_contractions_to_file(FILE *fout,int ncontr,int *op_sour,int *op_sink,complex *contr,int twall,
                                          const char *tag,double norm);

#endif
