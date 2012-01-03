#pragma once

///////////////// New types ///////////////////

typedef int coords[4];

typedef int intpair[2];
typedef uint32_t checksum[2];

typedef double complex[2];

typedef complex spin[4];
typedef complex color[3];

typedef spin colorspin[3];
typedef color spincolor[4];

typedef colorspin spincolorspin[4];
typedef spincolorspin colorspincolorspin[3];

typedef spin spinspin[4];
typedef spinspin colorspinspin[3];

typedef color su3[3];
typedef su3 quad_su3[4];

typedef colorspinspin su3spinspin[3];

typedef complex as2t[6];
typedef su3 as2t_su3[6];

//Random types
enum rnd_type{RND_ALL_PLUS_ONE,RND_ALL_MINUS_ONE,RND_UNIF,RND_Z2,RND_Z4,RND_GAUSS};
//Source type
enum source_type{POINT_SOURCE,UNDILUTED_SOURCE,COLOR_DILUTED_SOURCE,SPIN_DILUTED_SOURCE,SPINCOLOR_DILUTED_SOURCE};

//The structure for the random generator
typedef struct
{
  int idum;
  int idum2;
  int iv[ran2_ntab];
  int iy;
} rnd_gen;

//The structure for gamma matrix
typedef struct
{
  int pos[4];
  complex entr[4];
} dirac_matr;

//nissa vector
typedef struct
{
  int nel;
  int size_per_el;
  
  char tag[nissa_vect_string_length];
  char type[nissa_vect_string_length];
  
  char file[nissa_vect_string_length];
  int line;

  void *prev;
  void *next;
} nissa_vect;

//nissa file reader
typedef struct
{
  int open;
  int reading;
  int nbytes_per_site;
  char *buf;
  MPI_File *reader_file;
  LemonReader *lemon_reader;
} nissa_reader;
