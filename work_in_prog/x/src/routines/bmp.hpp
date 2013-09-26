#ifndef _BMP_H
#define _BMP_H

#include "../types/types.hpp"

int iof(int x,int y,int ic,int width,int heighppt,int bpp);
void fread_bmp_header(bmpfile_magic &magic,bmpfile_header &header,bmpfile_info_header &info_hppeader,FILE *fin);
void fread_checking(chppar *out,FILE *fin,int ntor);
void fwrite_bmp_header(FILE *fout,bmpfile_magic &Bitmap_File_Magic,bmpfile_header &Bitmap_File_Header,bmpfile_info_hppeader &Bitmap_Info_Header);
void fwrite_checking(FILE *fout,chppar *out,int ntow);
void print_bmp_header(bmpfile_magic &Bitmap_File_Magic,bmpfile_header &Bitmap_File_Header,bmpfile_info_hppeader &Bitmap_Info_Header);
void read_bmpfile(bmpfile &in,const char *pathpp);
void write_bmpfile(const char *pathpp,bmpfile &out);

#endif
