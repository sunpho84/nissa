#include "nissa.hpp"

using namespace nissa;

// int T,L;
// vector_remap_t *remapper;

// void index_from_Neo_to_lx(int &rank_out,int &iel_out,int iel_in,void *pars)
// {
//   int mu_ord[4]={0,3,2,1};
  
//   //decript from, cls order
//   int glb_site_sour=iel_in+rank*8*locVolh;
//   int shift_comp=glb_site_sour%2;
//   glb_site_sour/=2;
//   int mu=mu_ord[glb_site_sour%4];
//   glb_site_sour/=4;
//   glb_site_sour*=2;
//   if(glblx_parity(glb_site_sour)==0) glb_site_sour++;
  
//   //check
//   if(glb_site_sour>=glbVol) CRASH("%d>=%d impossible!",glb_site_sour,glbVol);
  
//   //get coords
//   const Coords g=glb_coord_of_glblx(glb_site_sour);
//   std::swap(g[1],g[3]); //only 1 and 3 must be switched
//   if(shift_comp) g[mu]=(g[mu]+glbSize[mu]-1)%glbSize[mu];
  
//   //get final results
//   get_loclx_and_rank_of_coord(&iel_out,&rank_out,g);
//   iel_out=iel_out*4+mu;
// }

// void conf_convert(char *outpath,char *inpath)
// {  
//   //open
//   FILE *fin=fopen(inpath,"r");
//   if(fin==NULL) CRASH("opening %s",inpath);
  
//   //read header
//   int nr;
//   coords temp;
//   nr=fread(temp,sizeof(coords),1,fin);
//   if(nr!=1) CRASH("did not success in reading");
//   double plaq;
//   nr=fread(&plaq,sizeof(double),1,fin);
//   if(nr!=1) CRASH("did not success in reading");
//   MPI_Barrier(MPI_COMM_WORLD);
  
//   //if needed convert the lattice size to check
//   if(!little_endian) change_endianness((uint32_t*)temp,(uint32_t*)temp,4);
//   if(temp[0]!=glbSize[0]||temp[1]!=glbSize[1]) CRASH("conf of size %dx%d, expecting %dx%d",temp[0],temp[1],glbSize[0],glbSize[1]);
  
//   //convert the plaquette
//   if(!little_endian) change_endianness(plaq);
//   plaq/=3;
  
//   //seek to the correct point
//   fseek(fin,24+rank*sizeof(quad_su3)*locVol,SEEK_SET);
//   MPI_Barrier(MPI_COMM_WORLD);

//   //read
//   char *buf=nissa_malloc("buf",locVol*sizeof(quad_su3),char);
//   nr=fread(buf,sizeof(quad_su3),locVol,fin);
//   if(nr!=locVol) CRASH("did not success in reading the conf");
//   MPI_Barrier(MPI_COMM_WORLD);
 
//   //if needed convert the endianess of the conf
//   if(!little_endian) change_endianness((double*)buf,(double*)buf,locVol*4*18);
  
//   //reorder
//   quad_su3 *conf=nissa_malloc("conf",locVol+bord_vol,quad_su3);
//   remapper->remap(conf,buf,sizeof(su3));

//   //compute the plaquette online
//   double plaq_comp=global_plaquette_lx_conf(conf);
//   MASTER_PRINTF("Plaquette computed: %16.16lg, read: %16.16lg\n",plaq,plaq_comp);
  
//   //write the conf
//   write_ildg_gauge_conf(outpath,conf,64);
  
//   //close
//   fclose(fin);
//   MPI_Barrier(MPI_COMM_WORLD);

//   nissa_free(buf);
//   nissa_free(conf);
// }

int main(int narg,char **arg)
{
  initNissa(narg,arg);
  CRASH(" ");
  
  // if(narg<2) CRASH("Use: %s input",arg[0]);
  
  // open_input(arg[1]);
  
  // // 1) Read information about the gauge conf
  
  // //Read the volume
  // read_str_int("L",&L);
  // read_str_int("T",&T);
  
  // //Init the MPI grid 
  // initGrid(T,L);
  
  // //init the remapper
  // remapper=new vector_remap_t(4*locVol,index_from_Neo_to_lx,NULL);
  
  // //read the number of gauge configurations
  // int N;
  // read_str_int("NGaugeConf",&N);
  
  // ///////////////////////////////////////////
  
  // for(int i=0;i<N;i++)
  // {
  //     char in[1024],out[1024];
  //     read_str(in,1024);
  //     read_str(out,1024);
      
  //     conf_convert(out,in);
  // }
  
  // close_input();
  
  // ///////////////////////////////////////////
  
  // delete remapper;
  
  closeNissa();
  
  return 0;
}
