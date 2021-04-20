#include "nissa.hpp"

using namespace nissa;

int T,L;
vector_remap_t *remapper;

void index_from_Neo_to_lx(Rank &rank_out,LocLxSite &iel_out,const LocLxSite& iel_in,void *pars)
{
  Coords<Dir> mu_ord;
  mu_ord(tDir)=tDir;
  mu_ord(xDir)=zDir;
  mu_ord(yDir)=yDir;
  mu_ord(zDir)=xDir;
  
  //decript from, cls order
  GlbLxSite glb_site_sour=iel_in()+rank*8*locVol()/2;
  const Parity shift_comp=(int)(glb_site_sour.nastyConvert()%2);
  glb_site_sour/=2;

  const Dir mu_sour=(int)(glb_site_sour()%4);
  const Dir mu=mu_ord(mu_sour);
  glb_site_sour/=4;
  glb_site_sour*=2;
  if(glblx_parity(glb_site_sour)==0) glb_site_sour++;
  
  //check
  if(glb_site_sour>=glbVol) crash("%d>=%d impossible!",glb_site_sour,glbVol);
  
  //get coords
  GlbCoords g;
  glb_coord_of_glblx(g,glb_site_sour);
  std::swap(g(xDir),g(zDir)); //only 1 and 3 must be switched
  if(shift_comp())
    g(mu)=(g(mu)+glbSize(mu)-1)%glbSize(mu);
  
  //get final results
  get_loclx_and_rank_of_coord(iel_out,rank_out,g);
  iel_out=iel_out*4+mu();
}

void conf_convert(char *outpath,char *inpath)
{  
  //open
  FILE *fin=fopen(inpath,"r");
  if(fin==NULL) crash("opening %s",inpath);
  
  //read header
  int nr;
  int temp[NDIM];
  nr=fread(temp,sizeof(int),NDIM,fin);
  if(nr!=1) crash("did not success in reading");
  double plaq;
  nr=fread(&plaq,sizeof(double),1,fin);
  if(nr!=1) crash("did not success in reading");
  MPI_Barrier(MPI_COMM_WORLD);
  
  //if needed convert the lattice size to check
  if(!little_endian) change_endianness((uint32_t*)temp,(uint32_t*)temp,4);
  if(temp[0]!=glbTimeSize or temp[1]!=glbSize(xDir))
    crash("conf of size %dx%d, expecting %dx%d",temp[0],temp[1],glbTimeSize(),glbSize(xDir)());
  
  //convert the plaquette
  if(!little_endian) change_endianness(plaq);
  plaq/=3;
  
  //seek to the correct point
  fseek(fin,24+rank*sizeof(quad_su3)*locVol(),SEEK_SET);
  MPI_Barrier(MPI_COMM_WORLD);

  //read
  char *buf=nissa_malloc("buf",locVol.nastyConvert()*sizeof(quad_su3),char);
  nr=fread(buf,sizeof(quad_su3),locVol(),fin);
  if(nr!=locVol) crash("did not success in reading the conf");
  MPI_Barrier(MPI_COMM_WORLD);
 
  //if needed convert the endianess of the conf
  if(!little_endian) change_endianness((double*)buf,(double*)buf,locVol.nastyConvert()*4*18);
  
  //reorder
  quad_su3 *conf=nissa_malloc("conf",locVolWithBord.nastyConvert(),quad_su3);
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
  remapper=new vector_remap_t(4*locVol(),index_from_Neo_to_lx,NULL);
  
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
