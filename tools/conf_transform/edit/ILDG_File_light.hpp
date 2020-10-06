#include <stdint.h>
#include <stdio.h>

typedef off_t ILDG_Offset;
typedef FILE* ILDG_File;

//ILDG header
struct ILDG_header
{
  uint32_t magic_no;
  uint16_t version;
  uint16_t mbme_flag;
  uint64_t data_length;
  char type[128];
};

void check_endianness();
ILDG_File ILDG_File_open(const char *path,const char *mode);
int ILDG_File_search_record(ILDG_header &header,ILDG_File &file,const char *record_name);
ILDG_Offset ILDG_File_get_position(ILDG_File &file);
void ILDG_File_read(void *data,ILDG_File &file,int nbytes_req);
void ILDG_File_close(ILDG_File &file);
void ILDG_File_set_position(ILDG_File &file,ILDG_Offset pos,int amode=SEEK_SET);
void ILDG_File_write(ILDG_File &file,void *data,int nbytes_req);
