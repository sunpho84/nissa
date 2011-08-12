#pragma once

///////////////// New types ///////////////////

typedef double complex[2];

typedef complex spin[4];
typedef complex color[3];

typedef spin colorspin[3];
typedef color redspincolor[2];
typedef color spincolor[4];

typedef spin spinspin[4];
typedef spinspin colorspinspin[3];

typedef color su3[3];
typedef su3 quad_su3[4];

typedef colorspinspin su3spinspin[3];

typedef complex as2t[6];
typedef su3 as2t_su3[6];

//The structure for gamma matrix
typedef struct
{
  int pos[4];
  complex entr[4];
} dirac_matr;

//Appretto vector
typedef struct
{
  int nel;
  int size_per_el;
  
  char tag[appretto_vect_string_length];
  char type[appretto_vect_string_length];
  
  char file[appretto_vect_string_length];
  int line;

  void *prev;
  void *next;
} appretto_vect;

//Appretto file reader
typedef struct
{
  int open;
  int reading;
  int nbytes_per_site;
  char *buf;
  MPI_File *reader_file;
  LemonReader *lemon_reader;
} appretto_reader;
