#pragma once

// ----------------------- double to float and vice-versa conversion ----------------------

void double_to_float_vec(float *out,double *in,int nel)
{for(int iel=0;iel<nel;iel++) out[iel]=(float)(in[iel]);}

void float_to_double_vec(double *out,float *in,int nel)
{for(int iel=0;iel<nel;iel++) out[iel]=(double)(in[iel]);}


// -------------------------- su3 matrix compression amd decompression ---------------------

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

//decompress 8 parameters to an su3 matrix
void su3c_to_su3f(su3f out,su3c in)
{
#include "gauge_decompression.c"
} 

//gpu version
__device__ void gpu_su3c_to_su3f(su3f out,su3c in)
{
#include "gauge_decompression.c"
} 


// ----------------------- full gaugeconf compression and decompression ----------------------

//compress an array of su3
void gaugeconf_compress(quad_su3c *out,quad_su3d *in,int size)
{for(int i=0;i<4*size;i++) su3d_to_su3c(((su3c*)out)[i],((su3d*)in)[i]);}

//reconstruction of a full gauge conf
__global__ void gpu_gaugeconf_reconstruct(quad_su3f *out,quad_su3c *in,int vol)
{
  int pos=threadIdx.x+blockDim.x*blockIdx.x;
  if(pos<vol*4)
    gpu_su3c_to_su3f(((su3f*)out)[pos],((su3c*)in)[pos]);
}

//Compress the configuration on cpu and send it to the gpu
void gaugeconf_compress_and_move(quad_su3c *gpu_comp_conf,quad_su3d *cpu_uncomp_conf,int vol)
{
  //allocate room for the compressed conf on cpu
  quad_su3c *cpu_comp_conf=cpu_allocate_quad_su3c(vol,"cpu_comp_conf");
  
  //compress the configuration and send it
  gaugeconf_compress(cpu_comp_conf,cpu_uncomp_conf,vol);
  memcpy_cpu_to_gpu(gpu_comp_conf,cpu_comp_conf,vol*sizeof(quad_su3c));
  
  //free cpu compressed conf
  free(cpu_comp_conf);
}

//Decompress the configuration on gpu and send it to the cpu
void gaugeconf_decompress_and_move(quad_su3f *cpu_uncomp_conf,quad_su3c *gpu_comp_conf,int vol,int blocksize)
{
  //allocate room for the decompressed conf on gpu
  quad_su3f *gpu_uncomp_conf=gpu_allocate_quad_su3f(vol,"gpu_uncomp_conf");
  
  //expand the configuration and move it
  int gridsize=(int)ceil(vol*4.0/blocksize);
  gpu_gaugeconf_reconstruct<<<gridsize,blocksize>>>(gpu_uncomp_conf,gpu_comp_conf,vol);
  memcpy_gpu_to_cpu(cpu_uncomp_conf,gpu_uncomp_conf,vol*sizeof(quad_su3f));

  //free the gpu decompressed conf
  gpu_free(gpu_uncomp_conf,"gpu_uncomp_conf");
}


// ------------------------------------ testing routines -----------------------------------

extern "C" void test_gaugeconf_compression(quad_su3d *conf,int vol)
{
  int blocksize=192;
  
  //Compress the configuration on cpu and send it to the gpu
  quad_su3c *gpu_comp_conf=gpu_allocate_quad_su3c(vol,"gpu_comp_conf");
  gaugeconf_compress_and_move(gpu_comp_conf,conf,vol);
  
  //Decompress the configuration on the gpu and get it back
  quad_su3f *cpu_uncomp_conf=cpu_allocate_quad_su3f(vol,"cpu_uncomp_conf");
  gaugeconf_decompress_and_move(cpu_uncomp_conf,gpu_comp_conf,vol,blocksize);
  
  //Now compare the decompressed conf with the original one
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
  
  //free the gpu compressed conf and cpu decompressed conf
  gpu_free(gpu_comp_conf,"gpu_comp_conf");
  free(cpu_uncomp_conf);
}

