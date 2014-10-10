#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>

#include "../../../../src/nissa.hpp"
using namespace std;

#include "../types/types.hpp"

int iof(int x,int y,int ic,int width,int height,int bpp)
{return (bpp>>3)*(x+y*width)+ic;}

void fread_checking(char *out,FILE *fin,int ntor)
{
  int nr=fread(out,1,ntor,fin);
  if(nr!=ntor)
    {
      perror("Error reading data");
      exit(1);
    }
}

void fwrite_checking(FILE *fout,char *out,int ntow)
{
  int nw=fwrite(out,1,ntow,fout);
  if(nw!=ntow)
    {
      perror("Error writing data");
      exit(1);
    }
}

void fread_bmp_header(bmpfile_magic &magic,bmpfile_header &header,bmpfile_info_header &info_header,FILE *fin)
{
  fread_checking((char*)(&magic),fin,sizeof(bmpfile_magic));
  fread_checking((char*)(&header),fin,sizeof(bmpfile_header));
  fread_checking((char*)(&info_header),fin,sizeof(bmpfile_info_header));
}

void print_bmp_header(bmpfile_magic &Bitmap_File_Magic,bmpfile_header &Bitmap_File_Header,bmpfile_info_header &Bitmap_Info_Header)
{
  printf("magic number %c%c\n",Bitmap_File_Magic.magic[0],Bitmap_File_Magic.magic[1]);
  printf("file_size %d\n",Bitmap_File_Header.filesz);
  printf("creator1 %d\n",Bitmap_File_Header.creator1);
  printf("creator2 %d\n",Bitmap_File_Header.creator2);
  printf("bmp_offset %d\n",Bitmap_File_Header.bmp_offset);
  
  printf("header_sz %d\n",Bitmap_Info_Header.header_sz);
  printf("width %d\n",Bitmap_Info_Header.width);
  printf("height %d\n",Bitmap_Info_Header.height);
  printf("nplanes %d\n",Bitmap_Info_Header.nplanes);
  printf("bitspp %d\n",Bitmap_Info_Header.bitspp);
  printf("compress_type %d\n",Bitmap_Info_Header.compress_type);
  printf("bmp_bytesz %d\n",Bitmap_Info_Header.bmp_bytesz);
  printf("hres %d\n",Bitmap_Info_Header.hres);
  printf("vres %d\n",Bitmap_Info_Header.vres);
  printf("ncolors %d\n",Bitmap_Info_Header.ncolors);
  printf("nimpcolors %d\n",Bitmap_Info_Header.nimpcolors);
}

void fwrite_bmp_header(FILE *fout,bmpfile_magic &Bitmap_File_Magic,bmpfile_header &Bitmap_File_Header,bmpfile_info_header &Bitmap_Info_Header)
{
  fwrite_checking(fout,(char*)(&Bitmap_File_Magic),sizeof(bmpfile_magic));
  fwrite_checking(fout,(char*)(&Bitmap_File_Header),sizeof(bmpfile_header));
  fwrite_checking(fout,(char*)(&Bitmap_Info_Header),sizeof(bmpfile_info_header));
}  

void read_bmpfile(bmpfile &in,const char *path)
{
  FILE *fin=open_file(path,"r");
  fread_bmp_header(in.magic,in.header,in.info_header,fin);
  print_bmp_header(in.magic,in.header,in.info_header);
  int size=in.info_header.bmp_bytesz;
  in.data=nissa_malloc("in",size,char);
  fread_checking(in.data,fin,size);
  fclose(fin);
}

void write_bmpfile(const char *path,bmpfile &out)
{
  FILE *fout=open_file(path,"w");
  fwrite_bmp_header(fout,out.magic,out.header,out.info_header);
  int size=out.info_header.bmp_bytesz;
  fwrite_checking(fout,out.data,size);
  fclose(fout);
}
