#pragma once

void double_to_float_vec(float *out,double *in,int nel)
{for(int iel=0;iel<nel;iel++) out[iel]=(float)(in[iel]);}

void float_to_double_vec(double *out,float *in,int nel)
{for(int iel=0;iel<nel;iel++) out[iel]=(double)(in[iel]);}

//compress an su3 matrix to 8 parameters
void su3d_to_su3c(su3c out,su3d in)
{
  out[0]=(float)(in[0][1][0]); //a2
  out[1]=(float)(in[0][1][1]);
  out[2]=(float)(in[0][2][0]); //a3
  out[3]=(float)(in[0][2][1]);
  out[4]=(float)(atan2(in[0][0][1],in[0][0][0])); //theta_a1
  out[5]=(float)(atan2(in[2][0][1],in[2][0][0])); //theta_a2
  out[6]=(float)(in[1][0][0]); //b1
  out[7]=(float)(in[1][0][1]);
}

//compress an array of su3
void gaugeconf_compress(quad_su3c *out,quad_su3d *in,int size)
{for(int i=0;i<4*size;i++) su3d_to_su3c(((su3c*)out)[i],((su3d*)in)[i]);}

//decompress 8 parameters to an su3 matrix
void su3c_to_su3f(su3f out,su3c in)
{
#include "gauge_decompression.c"
} 

//device version
__device__ void dev_su3c_to_su3f(su3f out,su3c in)
{
#include "gauge_decompression.c"
} 

//reconstruction of a full gauge conf
__global__ void gpu_gaugeconf_reconstruct(quad_su3f *out,quad_su3c *in,int vol)
{
  int pos=threadIdx.x+blockDim.x*blockIdx.x;
  if(pos<vol*4)
    dev_su3c_to_su3f(((su3f*)out)[pos],((su3c*)in)[pos]);
}

extern "C" void test_gaugeconf_compression(quad_su3d *conf,int vol)
{

  // -------- 1) Compress the configuration on cpu and send it to the gpu -----------

  //allocate room for the compressed conf on cpu
  int comp_size=vol*sizeof(quad_su3c);
  quad_su3c *cpu_comp_conf=(quad_su3c*)cpu_allocate_vector(comp_size,"cpu_comp_conf");
  
  //compress the configuration
  gaugeconf_compress(cpu_comp_conf,conf,vol);
  
  //allocate room for the compressed conf on the gpu and send it
  quad_su3c *gpu_comp_conf=(quad_su3c*)gpu_allocate_vector(comp_size,"gpu_comp_conf");
  memcpy_cpu_to_gpu(gpu_comp_conf,cpu_comp_conf,comp_size);
  
  //free cpu compressed conf
  free(cpu_comp_conf);
  
  
  // -------- 2) decompress the configuration on the gpu ------------
  
  //allocate room for the decompressed conf on gpu
  int uncomp_size=vol*sizeof(quad_su3f);
  quad_su3f *gpu_uncomp_conf=(quad_su3f*)gpu_allocate_vector(uncomp_size,"gpu_uncomp_conf");
  
  //expand the configuration
  int blocksize=192;
  int gridsize=(int)ceil(vol*4.0/blocksize);
  gpu_gaugeconf_reconstruct<<<gridsize,blocksize>>>(gpu_uncomp_conf,gpu_comp_conf,vol);

  //free gpu compressed conf
  gpu_free(gpu_comp_conf,"gpu_comp_conf");
  

  // -------- 3) copy the configuration back to cpu --------------
  
  //allocate room for the uncompressed conf on cpu and get it
  quad_su3f *cpu_uncomp_conf=(quad_su3f*)cpu_allocate_vector(uncomp_size,"cpu_uncomp_conf");
  memcpy_gpu_to_cpu(cpu_uncomp_conf,gpu_uncomp_conf,uncomp_size);
  
  //free the gpu decompressed conf
  gpu_free(gpu_uncomp_conf,"gpu_uncomp_conf");

  // -------- 4) now compare the decompressed conf with the original one ------- 

  for(int ivol=0;ivol<vol;ivol++)
    for(int idir=0;idir<4;idir++)
      {
	float err=0;
	
	for(int i=0;i<3;i++)
	  for(int j=0;j<3;j++)
	    for(int ri=0;ri<2;ri++)
	      {
		double diff=conf[ivol][idir][i][j][ri]-cpu_uncomp_conf[ivol][idir][i][j][ri];
		err+=diff*diff;
	      }
	
	printf("%d %d %g\n",ivol,idir,err);
      }
  
  //free the cpu decompressed conf
  free(cpu_uncomp_conf);
}

