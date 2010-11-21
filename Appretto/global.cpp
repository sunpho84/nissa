int T,L,loc_vol;
int loc_size[4],nproc_dir[4];
int proc_pos[4]={0,0,0,0};
int rank,rank_tot,cart_rank;

MPI_Comm cart_comm;

typedef double complex[2];

typedef complex spin[4];
typedef complex color[3];

typedef spin colorspin[3];
typedef color spincolor[4];
