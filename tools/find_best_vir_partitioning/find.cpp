#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <stdarg.h>
#include <cstdio>
#include <stdint.h>
#include <stdlib.h>
#include <map>
#include <vector>
#include <tr1/array>

#if defined BGQ && !defined BGQ_EMU
 #include <firmware/include/personality.h>
 #include <spi/include/kernel/location.h>
#endif

//crash reporting the expanded error message
void crash(const char *templ,...)
{
  //expand error message
  char mess[1024];
  va_list ap;
  va_start(ap,templ);
  vsprintf(mess,templ,ap);
  va_end(ap);
  
  fprintf(stderr,"ERROR \"%s\".\n",mess);
  exit(1);
}

//type to hold possible BGQ blocks
struct torus_grid_t
{
  int grid[5];
  bool is_torus[5];
  torus_grid_t(int a,int b,int c,int d)
  {
    grid[0]=a;
    grid[1]=b;
    grid[2]=c;
    grid[3]=d;
    grid[4]=2;
    
    is_torus[0]=a>2;
    is_torus[1]=b>2;
    is_torus[2]=c>2;
    is_torus[3]=d>2;
    is_torus[4]=1;
  }
  int get_N(){return grid[0]*grid[1]*grid[2]*grid[3]*grid[4];}
#if defined BGQ && !defined BGQ_EMU
  torus_grid_t()
  {
    Personality_t pers;
    Kernel_GetPersonality(&pers,sizeof(pers));
    
    //get size
    grid[0]=pers.Network_Config.Anodes;
    grid[1]=pers.Network_Config.Bnodes;
    grid[2]=pers.Network_Config.Cnodes;
    grid[3]=pers.Network_Config.Dnodes;
    grid[4]=pers.Network_Config.Enodes;
    for(int idir=0;idir<5;idir++) is_torus[idir]=ND_GET_TORUS(idir,pers.Network_Config.NetFlags);
  }
#else
  torus_grid_t() {}
#endif
};
torus_grid_t N64(2,2,4,2),N128(2,2,4,4),N256(4,2,4,4),N512(4,4,4,4),N1024(8,4,4,4),N2048(8,8,4,4);

struct coords5D_t: std::tr1::array<int,5> {};
struct coords4D_t: std::tr1::array<int,4> {};

struct mapping_t: std::map<int,coords5D_t>{};

//hold nranks per dir
struct rank_grid_t
{
  int grid[4];
  rank_grid_t(int a,int b,int c,int d)
  {
    grid[0]=a;
    grid[1]=b;
    grid[2]=c;
    grid[3]=d;
  }
  rank_grid_t() {*this=rank_grid_t(0,0,0,0);}
  rank_grid_t(coords4D_t &c){*this=rank_grid_t(c[0],c[1],c[2],c[3]);}
  int get_N(){return grid[0]*grid[1]*grid[2]*grid[3];}
};

//assigned dirs
struct assignement_t
{
  int assigned_torus_dirs[4][5];
  int nassigned_torus_dirs[4];
  void reset(){nassigned_torus_dirs[0]=nassigned_torus_dirs[1]=nassigned_torus_dirs[2]=nassigned_torus_dirs[3]=0;}
  assignement_t(){reset();}
  void set_from_way(int iway)
  {
    reset();
    for(int alpha=0;alpha<5;alpha++)
      {
	int mu=iway%4;
	assigned_torus_dirs[mu][nassigned_torus_dirs[mu]++]=alpha;
	iway>>=2;
      }
  }
  bool fit_the_rank(rank_grid_t &rank,torus_grid_t &torus)
  {
    bool out=true;
    for(int mu=0;mu<4;mu++)
      {
	int t=rank.grid[mu];
	//printf("    %d %d    ",mu,t);
	for(int i=0;i<nassigned_torus_dirs[mu];i++)
	  {
	    t/=torus.grid[assigned_torus_dirs[mu][i]];
	    //printf(" (%d %d %d)",assigned_torus_dirs[mu][i],torus.grid[assigned_torus_dirs[mu][i]],t);
	  }
	//printf("\n");
	out&=(t==1);
      }
    
    return out;
  }
  bool is_torus(torus_grid_t &torus)
  {
    bool out=true;
    for(int mu=0;mu<4;mu++) out&=(nassigned_torus_dirs[mu]!=1||torus.is_torus[assigned_torus_dirs[mu][0]]);
    return out;
  }
  mapping_t create_mapping(torus_grid_t &torus,rank_grid_t &rank)
  {
    mapping_t out;
    coords5D_t torus_coords;
    
    for(torus_coords[0]=0;torus_coords[0]<torus.grid[0];torus_coords[0]++)
      for(torus_coords[1]=0;torus_coords[1]<torus.grid[1];torus_coords[1]++)
	for(torus_coords[2]=0;torus_coords[2]<torus.grid[2];torus_coords[2]++)
	  for(torus_coords[3]=0;torus_coords[3]<torus.grid[3];torus_coords[3]++)
	    for(torus_coords[4]=0;torus_coords[4]<torus.grid[4];torus_coords[4]++)
	      {
		//find rank
		int rank_id=0;
		for(int mu=0;mu<4;mu++)
		  {
		    rank_id*=rank.grid[mu];
		    int rank_coords_mu=0;
		    bool sign=0;
		    for(int itorus_dir=0;itorus_dir<nassigned_torus_dirs[mu];itorus_dir++)
		      {
			//move according to sign direction
			int alpha=assigned_torus_dirs[mu][itorus_dir];
			int s=torus.grid[alpha],c=torus_coords[alpha];
			if(sign==0) c=s-1-c;
			rank_coords_mu=rank_coords_mu*s+c;
			
			//on odd rows, change direction
			if(torus_coords[alpha]%2==1) sign=!sign;
		      }
		    rank_id+=rank_coords_mu;
		  }
		
		//assign in the mapping
		out[rank_id]=torus_coords;
	      }
    
    return out;
  }
};

//factorize a number
int factorize(std::vector<int> &list,int N)
{
  int nfatt=0;
  int fatt=2;
  
  while(N>1)
    {
      int div=N/fatt;
      int res=N-div*fatt;
      if(res!=0) fatt++;
      else 
        {
          N=div;
          list.push_back(fatt);
          nfatt++;
        }
    }
  
  return nfatt;
}

struct valid_partition_lister_t
{
  //number of ranks in each direction for current partitioning
  coords4D_t R;
  
  int L[4],NR;
  int nfact,ncombo,icombo;
  std::vector<int> list_fact;
  bool factorize_rank;
  std::map<coords4D_t,int> old_partitionings;
  
  valid_partition_lister_t(coords4D_t sides,int NR) : NR(NR)
  {
    //copy L
    for(int mu=0;mu<4;mu++) L[mu]=sides[mu];
    
    //compute global and local volume
    uint64_t V=(uint64_t)L[0]*L[1]*L[2]*L[3];
    int LV=V/NR;
    
    //check that the global lattice is a multiple of the number of ranks
    //and that all directions can be made even, and one direction virtually parallelized
    if(V%NR) crash("global volume must be a multiple of ranks number");
    if((V/NR)%32!=0) crash("local size must be a multiple of 32");
    
    //factorize the local volume and the number of ranks
    std::vector<int> list_fact_LV,list_fact_NR;
    int nfact_LV=factorize(list_fact_LV,LV);
    int nfact_NR=factorize(list_fact_NR,NR);
    
    //if nfact_LV>=nfact_NR factorize the number of rank, otherwise the local volume
    //in the first case we find the best way to assign the ranks to different directions
    //in the second case we find how many sites per direction to assign to each rank
    factorize_rank=(nfact_LV>=nfact_NR);
    nfact=factorize_rank ? nfact_NR : nfact_LV;
    list_fact=factorize_rank ? list_fact_NR : list_fact_LV;
    
    //compute the number of combinations: this is given by 4^nfact
    ncombo=1;
    for(int ifact=0;ifact<nfact;ifact++) ncombo*=4;
    icombo=0;
  }
  bool find_next_valid_partitioning()
  {
    bool valid_partitioning=1;
    do
      {
	//reset number of ranks per dir
	R[0]=R[1]=R[2]=R[3]=1;
	
	//find the partioning corresponding to icombo
	int ifact=nfact-1;
	valid_partitioning=1;
	do
	  {
	    //find the direction: this is given by the ifact digit of icombo wrote in base 4
	    int mu=(icombo>>(2*ifact)) & 0x3;
	    
	    //if we are factorizing local lattice, rank factor is given by list_fact, otherwise L/list_fact
	    R[mu]*=list_fact[ifact];
            
	    //check that the total volume L is a multiple and it is larger than the number of proc
	    valid_partitioning=(L[mu]%R[mu]==0 && L[mu]>=R[mu]);
	    if(valid_partitioning) ifact--;
	  }
	while(valid_partitioning && ifact>=0);
	
	/*
	printf("temped partition: %d\n",icombo);
	for(int mu=0;mu<4;mu++) printf(" %d %d (%d)\n",mu,R[mu],L[mu]/R[mu]);
	*/
	
	bool one_direction_multiple_of_4=false;
	if(valid_partitioning)
	  for(int mu=0;mu<4;mu++)
	    {
	      //if we are factorizing reciprocal lattice, convert back to rank grid
	      if(!factorize_rank) R[mu]=L[mu]/R[mu];
	      //check that lattice size is even in all directions
	      valid_partitioning&=((L[mu]/R[mu])%2==0);
	      //take note if this dir is a multiple of 4
	      one_direction_multiple_of_4|=((L[mu]/R[mu])%4==0);
	    }
	
	//check that at least one direction is a multiple of 4
	valid_partitioning&=one_direction_multiple_of_4;
	
	//skip all remaining factorization using the same structure
	if(!valid_partitioning) icombo+=(ifact>1) ? 1<<(2*(ifact-1)) : 1;
	
	//printf("valid: %d, %d %d\n",valid_partitioning,icombo,ncombo);
	if(valid_partitioning)
	  {
	    valid_partitioning&=(old_partitionings.find(R)==old_partitionings.end());
	    if(valid_partitioning) old_partitionings[R]=icombo;
	    else icombo++;
	  }
      }
    while(!valid_partitioning && icombo<ncombo);
    
    return valid_partitioning;
  }
};

//check that the rank grid fit into the torus
assignement_t find_torus_assignement(torus_grid_t torus,rank_grid_t rank,int &iway)
{
  //check that grid and torus have the same size
  if(torus.get_N()!=rank.get_N()) crash("torus has %d elements, rank %d",torus.get_N(),rank.get_N());
  
  //there exist 4^5 way we can distribute torus into grid
  const int nways=4*4*4*4*4;
  
  //scan each of them
  assignement_t out;
  do
  {
    out.set_from_way(iway++);
    //for(int mu=0;mu<4;mu++)
    //{
    //printf("  ");
    //for(int i=0;i<out.nassigned_torus_dirs[mu];i++) printf("%d ",out.assigned_torus_dirs[mu][i]);
    //printf("\n");
    //}
    //printf("way %d fit the rank: %d, is torus: %d\n",iway-1,out.fit_the_rank(rank,torus),out.is_torus(torus));
  }
  while((!out.fit_the_rank(rank,torus) || !out.is_torus(torus)) && iway<nways);
  
  return out;
}

//return the torus containing the passed number of ranks
torus_grid_t find_torus(int &torus_size)
{
  if(torus_size>0)
    switch(torus_size)
      {
      case 64:  return N64;  break;
      case 128: return N128; break;
      case 256: return N256; break;
      case 512: return N512; break;
      case 1024: return N1024; break;
      case 2048: return N2048; break;
      default: crash("unknown partition to use for %d ranks",torus_size); return N64;break;
      }
  else
    {
      torus_grid_t out;
      printf("Describe the torus\n");
      torus_size=1;
      for(int mu=0;mu<5;mu++)
	{
	  printf("Insert grid_size[%d]: ",mu);
	  scanf("%d",&(out.grid[mu]));
	  int temp;
	  printf(" is_torus? ");
	  scanf("%d",&temp);
	  out.is_torus[mu]=temp;
	  torus_size*=out.grid[mu];
	}
	  printf("Fixing the total number of ranks to %d\n",torus_size);
      
      return out;
    }
}

void compute_border_border2_size(int &border,int &border2,coords4D_t L,torus_grid_t &torus,rank_grid_t &rank)
{
  //compute total border
  border=border2=0;
  for(int mu=0;mu<4;mu++)
    {
      //compute border in mu dir
      int border_mu=1;
      for(int nu=0;nu<4;nu++)
	if(mu!=nu)
	  border_mu*=L[nu]/rank.grid[nu];
      
      //summ to the toral border
      border+=border_mu;
      border2+=border_mu*border_mu;
    }
}

int main(int narg,char **arg)
{
  coords4D_t sides;
  torus_grid_t torus;
  
#if defined BGQ && !defined BGQ_EMU
  if(narg==2)
    {
      FILE *fin=fopen(arg[1],"r");
      if(fin==NULL) crash("error opening %s",arg[1]);
      
      size_t len=0;
      char *line=NULL;
      if(getline(&line,&len,fin)==-1) crash("reading");
      int L;
      sscanf(line,"L %d",&L);
      if(L>0)
	{
	  int T;
	  if(getline(&line,&len,fin)==-1) crash("reading");
	  sscanf(line,"T %d",&T);
	  sides[0]=T;
	  sides[1]=L;
	  sides[2]=L;
	  sides[3]=L;
	}
      else
	{
	  if(getline(&line,&len,fin)==-1) crash("reading");
	  if(sscanf(line,"LT %d",&(sides[0]))<=0) crash("reading LT");
	  if(getline(&line,&len,fin)==-1) crash("reading");
	  if(sscanf(line,"LX %d",&(sides[1]))<=0) crash("reading LX");
	  if(getline(&line,&len,fin)==-1) crash("reading");
	  if(sscanf(line,"LY %d",&(sides[2]))<=0) crash("reading LY");
	  if(getline(&line,&len,fin)==-1) crash("reading");
	  if(sscanf(line,"LZ %d",&(sides[3]))<=0) crash("reading LZ");
	}
      fclose(fin);
      free(line);
      printf("Finding partition for %d ranks\n",torus.get_N());
      for(int i=0;i<5;i++) printf(" torus[%d]: %d\n",i,torus.grid[i]);
    }
  else
#endif
    {
      printf("Insert L: ");
      int L;
      scanf("%d",&L);
      if(L>0)
	{
	  printf("Insert T: ");
	  int T;
	  scanf("%d",&T);
	  sides[0]=T;
	  sides[1]=L;
	  sides[2]=L;
	  sides[3]=L;
	}
      else
	{
	  printf("Insert LT: ");
	  scanf("%d",&sides[0]);
	  printf("Insert LX: ");
	  scanf("%d",&sides[1]);
	  printf("Insert LY: ");
	  scanf("%d",&sides[2]);
	  printf("Insert LZ: ");
	  scanf("%d",&sides[3]);
	}
      printf("Insert nranks: ");
      int nranks;
      scanf("%d",&nranks);
      torus=find_torus(nranks);
    }
  int nranks=torus.get_N();
  valid_partition_lister_t lister(sides,nranks);
  printf("Finding partition for T=%d x X=%d x Y=%d x Z=%d, %d ranks\n",sides[0],sides[1],sides[2],sides[3],nranks);
  for(int i=0;i<5;i++) printf(" size[%d]: %d, torus: %d\n",i,torus.grid[i],torus.is_torus[i]);
  
  //find the best rank assignement
  assignement_t assignement;
  rank_grid_t rank;
  int ntrue_valid=0;
  bool good;
  int bord=-1,bord2=-1;
  do
    {
      //advance along the list of good partitions up to finding a good one
      good=lister.find_next_valid_partitioning();
      if(good)
	{
	  printf("Valid partition combo %d: %d %d %d %d\n",lister.icombo,
		 lister.R[0],lister.R[1],lister.R[2],lister.R[3]);
	  printf("-----------------------------------------\n");
	  
	  //find a way to assign the torus
	  int iway=0;
	  rank_grid_t new_rank(lister.R);
	  assignement_t new_assignement=find_torus_assignement(torus,new_rank,iway);
	  
	  //check that we did not arrive to the maximal number of ways
	  const int nways=4*4*4*4*4;
	  if(iway!=nways)
	    {
	      printf(" Assignement obtained with way: %d\n",iway);
	      ntrue_valid++;
	      
	      //check if this is better than other
	      int new_bord,new_bord2;
	      compute_border_border2_size(new_bord,new_bord2,sides,torus,new_rank);
	      if(new_bord<bord||
		 (new_bord==bord&&
		  (new_bord2<bord2||
		   (new_bord2==bord2&&
		    (sides[0]/new_rank.grid[0]>sides[0]/rank.grid[0]))))||
		 bord==-1)
		{
		  printf(" Found new champion, %d(%d) beat %d(%d)\n",new_bord,new_bord2,bord,bord2);
		  bord=new_bord;
		  bord2=new_bord2; 
		  rank=new_rank;
		  assignement=new_assignement;
		}
	      else printf(" Defeated: %d(%d) loose %d(%d)\n",new_bord,new_bord2,bord,bord2);
	    }
	  else printf(" No valid assignement of torus direction to rank found\n");
	  
	  printf("\n");
	}
    }
  while(good);
  
  printf("==========================================================================\n");
  printf(" In total found %lu valid partitioning, true valid: %d\n",lister.old_partitionings.size(),ntrue_valid);
  printf("\n");
  
  //define local volume
  int LL[4]={sides[0]/rank.grid[0],sides[1]/rank.grid[1],sides[2]/rank.grid[2],sides[3]/rank.grid[3]};
  
  //print the ranks
  printf("Rank grid: %d x %d x %d x %d\n",rank.grid[0],rank.grid[1],rank.grid[2],rank.grid[3]);
  printf("Local volume: %d x %d x %d x %d\n",LL[0],LL[1],LL[2],LL[3]);
  printf("\n");
  
  //print the assignement of the torus direction to the rank
  printf(" Assignement:\n");
  for(int mu=0;mu<4;mu++)
    {
      printf(" [mu %d]: ",mu);
      for(int i=0;i<assignement.nassigned_torus_dirs[mu];i++)
	printf("%d ",assignement.assigned_torus_dirs[mu][i]);
      printf("\n");
    }
  printf("\n");
  
  //create a suffix for file names
  char path_suffix[50];
  sprintf(path_suffix,"%dranks_for_T%d_X%d_Y%d_Z%d",nranks,sides[0],sides[1],sides[2],sides[3]);
  
  //create name file for the mapping and open it
  char path[50];
  sprintf(path,"mapping_file_%s",path_suffix);
  FILE *fout=fopen(path,"w");
  if(fout==NULL) crash("opening file %s",path);
  
  //print the mapping
  mapping_t mapping=assignement.create_mapping(torus,rank);
  for(int irank=0;irank<nranks;irank++)
    {
      for(int alpha=0;alpha<5;alpha++)
	fprintf(fout," %d",mapping[irank][alpha]);
      fprintf(fout," 0\n");
    }
  fclose(fout);
  
  //find the direction to use to virtual parallelize
  int v=-1;
  for(int mu=0;mu<4;mu++) if(LL[mu]%4==0&&(v==-1||LL[mu]>LL[v])) v=mu;
  
  //create name file for nissa_config
  sprintf(path,"nissa_config_%s",path_suffix);
  fout=fopen(path,"w");
  if(fout==NULL) crash("opening file %s",path);
  
  //create nissa_config file
  fprintf(fout,"vnode_paral_dir %d\n",v);
  const char tag[]="txyz";
  for(int mu=0;mu<4;mu++)
    fprintf(fout,"set_%c_nranks %d\n",tag[mu],rank.grid[mu]);
  fprintf(fout,"\n");
  fclose(fout);
  
  return 0;
}
