#include "nissa.hpp"

using namespace nissa;

int T,L;
vector_remap_t *remapper;

void index_from_Neo_to_lx(int &rank_out,int &iel_out,int iel_in,void *pars)
{
  //decript fro cls order
  int glb_site_sour=iel_in+rank*8*loc_volh;
  int shift_comp=glb_site_sour%2;
  glb_site_sour/=2;
  int mu=glb_site_sour%4;
  glb_site_sour/=4;
  glb_site_sour*=2;
  
  //check
  if(glb_site_sour>=glb_vol) crash("%d>=%d impossible!",glb_site_sour,glb_vol);
  
  //get coords
  coords g;
  glb_coord_of_glblx(g,glb_site_sour);
  std::swap(g[1],g[3]); //only 1 and 3 must be switched
  if(shift_comp) g[mu]=(g[mu]+glb_size[mu]-1)%glb_size[mu];
  
  //get final results
  get_loclx_and_rank_of_coord(&iel_out,&rank_out,g);
  iel_out=iel_out*4+mu;
}

void conf_convert(char *outpath,char *inpath)
{  
  //read the header
  FILE *fin=fopen(inpath,"r");
  
  //read header
  int nr;
  coords temp;
  nr=fread(temp,sizeof(coords),1,fin);
  if(nr!=1) crash("did not success in reading");
  double plaq;
  nr=fread(&plaq,sizeof(double),1,fin);
  if(nr!=1) crash("did not success in reading");
  MPI_Barrier(MPI_COMM_WORLD);
  
  //if needed convert the lattice size to check
  if(!little_endian) change_endianness((uint32_t*)temp,(uint32_t*)temp,4);
  if(temp[0]!=glb_size[0]||temp[1]!=glb_size[1]) crash("conf of size %dx%d, expecting %dx%d",temp[0],temp[1],glb_size[0],glb_size[1]);
  
  //convert the plaquette
  if(!little_endian) change_endianness(plaq);
  
  //seek to the correct point
  fseek(fin,24+rank*sizeof(quad_su3)*loc_vol,SEEK_SET);
  MPI_Barrier(MPI_COMM_WORLD);

  //read
  char *buf=nissa_malloc("buf",loc_vol*sizeof(quad_su3),char);
  nr=fread(buf,sizeof(quad_su3),loc_vol,fin);
  if(nr!=loc_vol) crash("did not success in reading the conf");
  MPI_Barrier(MPI_COMM_WORLD);
 
  //if needed convert the endianess of the conf
  if(!little_endian) change_endianness((double*)buf,(double*)buf,loc_vol*4*18);
  
  //reorder
  quad_su3 *conf=nissa_malloc("conf",loc_vol+bord_vol,quad_su3);
  remapper->remap(conf,buf,sizeof(su3));

  //compute the plaquette online
  double plaq_comp=global_plaquette_lx_conf(conf);
  master_printf("Plaquette computed: %16.16lg, read: %16.16lg\n",plaq,plaq_comp);
  
  //write the conf
  write_ildg_gauge_conf(outpath,conf,64);
  
  //close
  fclose(fin);
  MPI_Barrier(MPI_COMM_WORLD);

  nissa_free(buf);
  nissa_free(conf);
}

void in_main(int narg,char **arg)
{
  if(narg<2) crash("Use: %s input",arg[0]);
  
  open_input(arg[1]);
  
  // 1) Read information about the gauge conf
  
  //Read the volume
  read_str_int("L",&L);
  read_str_int("T",&T);
  
  //Init the MPI grid 
  init_grid(T,L);
  
  //init the remapper
  remapper=new vector_remap_t(4*loc_vol,index_from_Neo_to_lx,NULL);
  
  //read the number of gauge configurations
  int N;
  read_str_int("NGaugeConf",&N);
  
  ///////////////////////////////////////////
  
  for(int i=0;i<N;i++)
  {
      char in[1024],out[1024];
      read_str(in,1024);
      read_str(out,1024);
      
      conf_convert(out,in);
  }
  
  close_input();
  
  ///////////////////////////////////////////
  
  delete remapper;
}

int main(int narg,char **arg)
{
  init_nissa_threaded(narg,arg,in_main);
  close_nissa();
  
  return 0;
}
