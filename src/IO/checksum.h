#ifndef _CHECKSUM_H
#define _CHECKSUM_H

uint32_t ildg_crc32(uint32_t crc,const unsigned char *buf,size_t len);
uint32_t ildg_crc32_fix_endianess(uint32_t crc,const unsigned char *buf,size_t len);
void checksum_compute_ildg_data(uint32_t *check,void *data,size_t bps);
void checksum_compute_nissa_data(uint32_t *check,void *data,size_t bps);

#endif
